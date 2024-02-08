"""
from: http://stackoverflow.com/questions/1112343/how-do-i-capture-sigint-in-python
"""
import signal


class SignalHandler(object):
    def __init__(self, sigs=[signal.SIGINT, signal.SIGTERM]):
        self.signals = sigs

    def __enter__(self):

        self.interrupted = False
        self.released = False
        self.sig_used = None

        self.original_handlers = {}
        for sig in self.signals:
            self.original_handlers[sig] = signal.getsignal(sig)

        def handler(signum, frame):
            self.release()
            self.interrupted = True
            self.sig_used = signum

        for sig in self.signals:
            signal.signal(sig, handler)

        return self

    def __exit__(self, type, value, tb):
        self.release()

    def release(self):

        if self.released:
            return False

        for sig, handler in self.original_handlers.items():
            signal.signal(sig, handler)

        self.released = True

        return True

if __name__ == '__main__':

    import unittest
    import time

    class SignalHandlerTestCase(unittest.TestCase):

        def test_simple(self):
            with SignalHandler() as h:
                while True:
                    print("h")
                    time.sleep(1)
                    print("e")
                    time.sleep(1)
                    print("l")
                    time.sleep(1)
                    print("l")
                    time.sleep(1)
                    print("o")
                    time.sleep(1)
                    if h.interrupted:
                        print("interrupted!")
                        break


    unittest.main()
