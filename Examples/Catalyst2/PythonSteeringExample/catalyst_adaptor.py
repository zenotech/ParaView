import os
import glob
import sys
import catalyst
import catalyst_conduit
from mpi4py import MPI


def initialize():
    # initialize ParaView Catalyst V2.
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    node = catalyst_conduit.Node()
    count = 0
    for i in sys.argv[1:]:
        node['catalyst/scripts/script'+str(count)] = i
        count = count + 1
        if rank == 0:
            print('Using Catalyst script', i)
    node['catalyst_load/implementation'] = 'paraview'
    # if we know the location of libcatalyst-paraview.so and want to pass that along via
    # Python we would do that like:
    # node['catalyst_load/search_paths/paraview'] = </full/path/to/libcatalyst-paraview.so>

    catalyst.initialize(node)

def finalize():
    # finalize ParaView Catalyst V2.
    node = catalyst_conduit.Node()
    catalyst.finalize(node)

def results():
    # get results from catalyst script
    node = catalyst_conduit.Node()
    catalyst.results(node)
    return node

def execute(time, timeStep, grid, attributes, results = None):
    # do the actual in situ analysis and visualization.
    node = catalyst_conduit.Node()

    node['catalyst/state/timestep'] = timeStep
    node['catalyst/state/time'] = time

    # the Catalyst channel is "input"
    node['catalyst/channels/input/type'] = 'mesh'

    mesh = node['catalyst/channels/input/data']

    mesh['coordsets/coords/type'] = 'uniform'

    mesh['coordsets/coords/dims/i'] = grid.XEndPoint - grid.XStartPoint + 1
    mesh['coordsets/coords/dims/j'] = grid.NumberOfYPoints
    mesh['coordsets/coords/dims/k'] = grid.NumberOfZPoints

    mesh['coordsets/coords/origin/x'] = grid.XStartPoint * grid.Spacing[0]
    mesh['coordsets/coords/origin/y'] = 0.0
    mesh['coordsets/coords/origin/z'] = 0.0

    mesh['coordsets/coords/spacing/dx'] = grid.Spacing[0]
    mesh['coordsets/coords/spacing/dy'] = grid.Spacing[1]
    mesh['coordsets/coords/spacing/dz'] = grid.Spacing[2]

    mesh['topologies/mesh/type'] = 'uniform'
    mesh['topologies/mesh/coordset'] = 'coords'

    fields = mesh['fields']

    use_velocity = results is None or results['input/request/velocity'] == True
    use_pressure = results is None or results['input/request/pressure'] == True


    # velocity is point-data.
    if use_velocity:
      # velocity is stored in interlaced form.
      fields['velocity/association'] = 'vertex'
      fields['velocity/topology'] = 'mesh'
      fields['velocity/volume_dependent'] = 'false'
      # flatten creates a deep copy of the array so we need to deep copy
      # this data into the Conduit Node.
      fields['velocity/values/x'] = attributes.Velocity.flatten()[0::3]
      fields['velocity/values/y'] = attributes.Velocity.flatten()[1::3]
      fields['velocity/values/z'] = attributes.Velocity.flatten()[2::3]

    if use_pressure:
      # pressure is cell-data.
      fields['pressure/association'] = 'element'
      fields['pressure/topology'] = 'mesh'
      fields['pressure/volume_dependent'] = 'false'
      # pressure is zero-copied into Catalyst because the simulation
      # stores the data the same way that Conduit expects it and it is scalar.
      fields['pressure/values'].set_external(attributes.Pressure)


    catalyst.execute(node)

    return
