def stripargument(function):
    """Allow calling a one-argument function with two arguments."""
    
    if function.func_code.co_argcount == 1:
        return (lambda x,y: function(x))
    else:
        return function

def select(parameters,selected):
    if not isinstance(selected,(tuple,list)):
        selected = [selected]
    
    if isinstance(parameters,(tuple,list)):
        return [par for par in parameters if par in selected]
    elif isinstance(parameters,dict):
        return dict((key,parameters[key]) for key in selected)
    else:
        return TypeError

