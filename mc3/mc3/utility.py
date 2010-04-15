def stripargument(function):
    if function.func_code.co_argcount == 1:
        return (lambda x,y: function(x))
    else:
        return function
