# ****************************************************************************
#
# Copyright (C) 2019-2023, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
"""
import argparse
import sys
import json

from http.server import BaseHTTPRequestHandler, HTTPServer
from io import BytesIO



class RequestHandler(BaseHTTPRequestHandler):
    @cached_property
    def url(self):
        return urlparse(self.path)

    @cached_property
    def query_data(self):
        return dict(parse_qsl(self.url.query))

    @cached_property
    def post_data(self):
        content_length = int(self.headers.get("Content-Length", 0))
        return self.rfile.read(content_length)

    @cached_property
    def form_data(self):
        return dict(parse_qsl(self.post_data.decode("utf-8")))

    @cached_property
    def cookies(self):
        return SimpleCookie(self.headers.get("Cookie"))

    def do_GET(self):

        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        response = BytesIO()
        response.write(json.dumps({'message': 'GET request received'}).encode('utf-8'))
        self.wfile.write(response.getvalue())

    def do_POST(self):

        content_length = int(self.headers['Content-Length'])
        body = self.rfile.read(content_length)
        data = json.loads(body.decode('utf-8'))

        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        response = BytesIO()
        response.write(json.dumps({'message': 'POST request received', 'data_received': data}).encode('utf-8'))
        self.wfile.write(response.getvalue())

def run_alert_server(address='127.0.0.1', port=8080):
    """
    """
    httpd = HTTPServer((address, port), RequestHandler)

    try:
        httpd.serve_forever()
        print(f'Alert server started on port {port}')

    except KeyboardInterrupt:
        httpd.server_close()
        print(f'Alert server closed by user')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Run the alert server")

    parser.add_argument("-p", "--port", type=int,
                        default=8080,
                        help="Port number (default: 8080)")

    parser.add_argument("-a", "--address", type=str,
                        default="127.0.0.1",
                        help="Server address (default: 127.0.0.1)")
    
    args = parser.parse_args()

    run_alert_server(ip=args.address, port=args.port)

    
