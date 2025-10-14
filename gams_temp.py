from gamspy import Container, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

m = Container()

i = Set(Container=m, name='comps', description="Involved chemical components")
j = Set(Container=m, name='streams', description="Flow streams in process")

v1 = Parameter(
     Container=m,
     name='v1',
     domain=i,  # over components
     description="Stochiometric coefficient \
                  of component i in haber reactor"
)

F = Variable(
     Container=m,
     name="F",
     domain=[i, j],  # per comp & stream
     type="Positive",
     description="Molar flowrate of component i in stream j"
)