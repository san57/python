def execute(self, **kwargs):
    print(__file__)
    import code

    code.interact(local=dict(locals(), **globals()))
