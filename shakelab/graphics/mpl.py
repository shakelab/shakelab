import matplotlib.pyplot as plt
from shakelab.signals.base import Record, Stream, StreamCollection


def streamplot(obj):
    """
    Plot seismic traces from Record, Stream, or StreamCollection.

    Parameters:
    obj: Record, Stream, or StreamCollection instance.
    """

    if isinstance(obj, Record):
        fig, ax = plt.subplots()
        ax.plot(obj.taxis, obj.data)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude')
        ax.set_title('Seismic Trace')
        ax.legend()
        plt.show()

    elif isinstance(obj, Stream):
        fig, ax = plt.subplots()
        for record in obj.record:
            ax.plot(record.taxis, record.dat)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude')
        ax.set_title('Seismic Stream')
        ax.legend()
        plt.show()

    elif isinstance(obj, StreamCollection):
        nstreams = len(obj.stream)
        fig, axs = plt.subplots(nstreams, 1, sharex=True, figsize=(8, 2 * nstreams))
        if nstreams == 1:
            axs = [axs]

        for ax, stream in zip(axs, obj.stream):
            for record in stream.record:
                ax.plot(record.taxis, record.data)
            ax.set_ylabel('Amplitude')
            ax.legend()

        axs[-1].set_xlabel('Time (s)')
        fig.suptitle('Seismic StreamCollection')
        plt.tight_layout()
        plt.show()

    else:
        raise TypeError("Input must be Record, Stream, or StreamCollection.")