# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
HTTP requests to FDSN web services.

The requests use the standard library only (:mod:`urllib.request`):

- GET with the parameters in the URL, POST with a text body (selection
  lines, see :mod:`~shakelab.signals.fdsnws.selection`); a POST
  redirected with status 307 or 308 is sent again to the new address;
- the replies are requested compressed (gzip) and decompressed while
  they are read, also when written to a file block by block;
- certificates are verified with those named by the environment
  variables ``SSL_CERT_FILE`` or ``SSL_CERT_DIR``, if set, otherwise
  with the bundle of :mod:`certifi`, if installed, otherwise with
  those of the system; the verification is never disabled;
- timeouts, connection failures and the statuses 429, 500, 502, 503
  and 504 are retried with growing waits; other statuses raise
  :class:`FdsnError` with the message of the server; 204 means no data.

The parameter values follow the FDSN conventions (Commonalities 1.2):
times as ``YYYY-MM-DDThh:mm:ss.ssssss`` in UTC, real numbers in decimal
notation, booleans as ``true``/``false``.
"""

from __future__ import annotations

import gzip
import os
import socket
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from decimal import Decimal
from numbers import Integral, Real
from typing import Any, Optional

from shakelab.libutils.timeN import Date

__all__ = [
    "DEFAULT_RETRIES",
    "DEFAULT_TIMEOUT",
    "FdsnError",
    "FdsnWarning",
    "format_value",
    "request",
    "user_agent",
]

DEFAULT_TIMEOUT = 120.0
DEFAULT_RETRIES = 3

# Statuses worth another attempt, and the first wait (seconds; doubled
# at each attempt, at most RETRY_MAX_WAIT).
RETRY_STATUSES = (429, 500, 502, 503, 504)
RETRY_WAIT = 1.0
RETRY_MAX_WAIT = 30.0

_CHUNK = 1 << 16
_CONTEXT: Optional[ssl.SSLContext] = None


# ----------------------------------------------------------------------------
# Errors
# ----------------------------------------------------------------------------


class FdsnError(RuntimeError):
    """
    Error of an FDSN web service request.

    Parameters
    ----------
    message : str
        Description of the error.
    status : int, optional
        HTTP status of the reply, if any.
    url : str, optional
        Address of the request.
    detail : str, optional
        Text of the reply of the server (its error message).
    """

    def __init__(
        self,
        message: str,
        status: Optional[int] = None,
        url: Optional[str] = None,
        detail: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.status = status
        self.url = url
        self.detail = detail


class FdsnWarning(UserWarning):
    """Warning about the data returned by an FDSN web service."""


# ----------------------------------------------------------------------------
# Requests
# ----------------------------------------------------------------------------


def request(
    url: str,
    params: Optional[dict] = None,
    body: Optional[str] = None,
    target: Any = None,
    timeout: Real = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
) -> Any:
    """
    Send a request to a web service and return the reply.

    Parameters
    ----------
    url : str
        Address of the method (e.g. ``.../fdsnws/station/1/query``).
    params : dict, optional
        Parameters of a GET request, formatted with
        :func:`format_value`; None values are left out.
    body : str, optional
        Body of a POST request; if given, ``params`` must be None.
    target : str or path-like, optional
        File where the reply is written, through a temporary file
        (``<target>.part``) renamed at the end, so that no partial file
        is left; an existing file is replaced.
    timeout : float, optional
        Seconds without data before a request is abandoned.
    retries : int, optional
        Attempts after the first one, for the failures listed in the
        module notes.

    Returns
    -------
    bytes, int or None
        The reply, or the number of bytes written to ``target``; None
        if the service has no data (status 204).

    Raises
    ------
    FdsnError
        If the service answers with an error, or if it cannot be
        reached after the retries.
    """
    if params is not None and body is not None:
        raise ValueError("a request has either parameters or a body")
    if params:
        query = urllib.parse.urlencode(
            {
                key: format_value(value)
                for key, value in params.items()
                if value is not None
            },
            safe="*?,:",
        )
        url = f"{url}?{query}"
    data = None if body is None else body.encode("utf-8")

    attempt = 0
    while True:
        try:
            return _send(url, data, target, timeout)
        except _Retry as error:
            if attempt >= retries:
                raise error.error from None
            time.sleep(min(error.wait(attempt), RETRY_MAX_WAIT))
            attempt += 1


def user_agent() -> str:
    """
    Return the User-Agent of the requests.

    Returns
    -------
    str
        ``ShakeLab/<version> (https://shakelab.org)``.
    """
    from shakelab import __version__

    return f"ShakeLab/{__version__} (https://shakelab.org)"


def format_value(value: Any) -> str:
    """
    Format a parameter value as the FDSN web services expect.

    Parameters
    ----------
    value : Date, bool, int, float or str
        Value to format.

    Returns
    -------
    str
        Times in UTC with microseconds and no zone designator, booleans
        as ``true`` or ``false``, real numbers in decimal notation
        (never exponential), strings unchanged.

    Raises
    ------
    TypeError
        If the value is of another type.
    ValueError
        If a real number is not finite.

    Examples
    --------
    >>> format_value(1e-05)
    '0.00001'
    >>> format_value(True)
    'true'
    """
    if isinstance(value, Date):
        return value.to_iso8601().rstrip("Z")
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, Integral):
        return str(int(value))
    if isinstance(value, Real):
        value = float(value)
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError("parameter values must be finite")
        text = format(Decimal(repr(value)), "f")
        if "." in text:
            text = text.rstrip("0").rstrip(".")
        return text or "0"
    if isinstance(value, str):
        return value
    raise TypeError(f"unsupported parameter value: {value!r}")


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


class _Retry(Exception):
    """A failure worth another attempt, with the error to raise last."""

    def __init__(self, error: FdsnError, retry_after: Optional[float]):
        super().__init__(str(error))
        self.error = error
        self.retry_after = retry_after

    def wait(self, attempt: int) -> float:
        if self.retry_after is not None:
            return self.retry_after
        return RETRY_WAIT * 2 ** attempt


class _RedirectHandler(urllib.request.HTTPRedirectHandler):
    """Redirects that keep the method and body of 307 and 308 replies."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code in (307, 308) and req.get_method() == "POST":
            return urllib.request.Request(
                newurl,
                data=req.data,
                headers=dict(req.header_items()),
                method="POST",
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _context() -> ssl.SSLContext:
    """
    TLS context (cached): the certificates named by ``SSL_CERT_FILE``
    or ``SSL_CERT_DIR`` if set (e.g. behind a proxy), otherwise the
    certifi bundle if installed, otherwise those of the system.
    """
    global _CONTEXT
    if _CONTEXT is None:
        certifi = None
        if not (
            os.environ.get("SSL_CERT_FILE") or os.environ.get("SSL_CERT_DIR")
        ):
            try:
                import certifi
            except ImportError:
                certifi = None
        if certifi is None:
            _CONTEXT = ssl.create_default_context()
        else:
            _CONTEXT = ssl.create_default_context(cafile=certifi.where())
    return _CONTEXT


def _opener() -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(
        _RedirectHandler,
        urllib.request.HTTPSHandler(context=_context()),
    )


def _send(url: str, data: Optional[bytes], target: Any, timeout: Real):
    """One attempt; raise _Retry for the failures worth retrying."""
    headers = {"User-Agent": user_agent(), "Accept-Encoding": "gzip"}
    if data is not None:
        headers["Content-Type"] = "text/plain"
    req = urllib.request.Request(url, data=data, headers=headers)

    try:
        reply = _opener().open(req, timeout=timeout)
    except urllib.error.HTTPError as error:
        status = error.code
        try:
            detail = _error_text(error)
            retry_after = _retry_after(error.headers)
        finally:
            # The error holds the open reply: close it now, not when
            # it is garbage collected.
            error.close()
        message = f"HTTP {status} from {url}"
        if detail:
            message = f"{message}: {detail.splitlines()[0][:200]}"
        if status == 413:
            message += " (request too large: split it, e.g. with chunk=)"
        failure = FdsnError(message, status=status, url=url, detail=detail)
        if status in RETRY_STATUSES:
            raise _Retry(failure, retry_after) from None
        raise failure from None
    except urllib.error.URLError as error:
        reason = error.reason
        if isinstance(reason, ssl.SSLCertVerificationError):
            raise FdsnError(
                f"certificate verification failed for {url}: {reason}. "
                "Install the certifi package, or the certificates of "
                "your Python installation (on macOS, with the Python of "
                "python.org: 'Install Certificates.command').",
                url=url,
            ) from None
        failure = FdsnError(f"cannot reach {url}: {reason}", url=url)
        raise _Retry(failure, None) from None
    except (TimeoutError, socket.timeout) as error:
        failure = FdsnError(f"timeout from {url}: {error}", url=url)
        raise _Retry(failure, None) from None

    with reply:
        if reply.status == 204:
            return None
        gzipped = reply.headers.get("Content-Encoding", "") == "gzip"
        try:
            if target is None:
                content = reply.read()
                return gzip.decompress(content) if gzipped else content
            return _write(reply, target, gzipped)
        except (TimeoutError, socket.timeout, ConnectionError) as error:
            failure = FdsnError(
                f"transfer from {url} interrupted: {error}", url=url
            )
            raise _Retry(failure, None) from None


def _write(reply: Any, target: Any, gzipped: bool) -> int:
    """Write a reply to a file through a temporary file."""
    path = os.fspath(target)
    partial = f"{path}.part"
    try:
        with open(partial, "wb") as fobj:
            size = _copy(reply, fobj, gzipped)
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.remove(partial)
    return size


def _copy(reply: Any, fobj: Any, gzipped: bool) -> int:
    """Copy a reply to a file block by block, decompressing it."""
    inflater = zlib.decompressobj(16 + zlib.MAX_WBITS) if gzipped else None
    size = 0
    while True:
        block = reply.read(_CHUNK)
        if not block:
            break
        if inflater is not None:
            block = inflater.decompress(block)
        fobj.write(block)
        size += len(block)
    if inflater is not None:
        block = inflater.flush()
        fobj.write(block)
        size += len(block)
    return size


def _error_text(error: urllib.error.HTTPError) -> str:
    """Text of an error reply, decompressed, without blank lines."""
    try:
        content = error.read()
    except OSError:
        return ""
    if error.headers.get("Content-Encoding", "") == "gzip":
        try:
            content = gzip.decompress(content)
        except OSError:
            return ""
    text = content.decode("utf-8", errors="replace")
    return "\n".join(line for line in text.splitlines() if line.strip())


def _retry_after(headers: Any) -> Optional[float]:
    """Seconds of a Retry-After header, if given in seconds."""
    value = headers.get("Retry-After") if headers is not None else None
    try:
        return max(0.0, float(value))
    except (TypeError, ValueError):
        return None
