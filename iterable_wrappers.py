class Simple:
    def __init__(self, iterable, *args, **kwargs):
        self.iterable = iterable
        if len(args) > 0:
            print(args[0] + "...")
        elif "desc" in kwargs:
            print(kwargs["desc"] + "...")

    def __iter__(self):
        return iter(self.iterable)


class Quiet:
    def __init__(self, iterable, *_, **__):
        self.iterable = iterable

    def __iter__(self):
        return iter(self.iterable)