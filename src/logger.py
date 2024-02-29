class Logger:
    def __init__(self):
        self.msgs = []

    def _log(self, msg, end="\n"):
        self.msgs.append(msg + end)
        print(msg)

    def _save_log(self, path):
        with open(path, 'w') as fd:
            for msg in self.msgs:
                fd.write(msg)
