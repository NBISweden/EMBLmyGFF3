#!/usr/bin/env python3
"""
ConciseStreamHandler is a subclass of logging.StreamHandler which will maintain
a list of previously emitted messages in order to have a more readable log.
"""

import logging

class ConciseStreamHandler(logging.StreamHandler):
    """
    ConciseStreamHandler is a subclass of logging.StreamHandler which will
    maintain a list of previously emitted messages in order to have a more
    readable log.
    """
    def __init__(self):
        super().__init__()
        self.log = {}

    def emit(self, record):
        try:
            msg = self.format(record)

            if record.message in self.log:
                self.log[record.message] += 1
            else:
                self.log[record.message] = 1
                self.stream.write(msg)
                self.stream.write(self.terminator)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)
