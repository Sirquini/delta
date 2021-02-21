
def dynamic_class(normal_class):
    class Wrapper(normal_class):
        @classmethod
        def set_method(cls, method):
            assert hasattr(method, '__call__'), f'Not callable method: {method}'
            setattr(cls, method.__name__, method)

        @classmethod
        def set_classmethod(cls, method):
            assert hasattr(method, '__call__'), f'Not callable method: {method}'
            setattr(cls, method.__name__, classmethod(method))

        @classmethod
        def set_staticmethod(cls, method):
            assert hasattr(method, '__call__'), f'Not callable method: {method}'
            setattr(cls, method.__name__, staticmethod(method))

        @classmethod
        def set_property(cls, method):
            assert hasattr(method, '__call__'), f'Not callable method: {method}'
            setattr(cls, method.__name__, property(method))

    Wrapper.__name__=f'dynamic_class({normal_class.__name__})'
    Wrapper.__doc__=(
        f'    Dynamic class of {normal_class.__name__} with three decorators:\n'
        f'        @{normal_class.__name__}.set_method\n'
        f'        @{normal_class.__name__}.set_classmethod\n'
        f'        @{normal_class.__name__}.set_staticmethod\n'
        f'        @{normal_class.__name__}.set_property\n'
        f'    \n'
        f'    Documentation of {normal_class.__name__}:\n'
        f'{normal_class.__doc__}'
    )
    return Wrapper



