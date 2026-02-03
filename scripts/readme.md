This is a helper script to build topology (martini3) for mixed fibrils. You need to have run pdb (triple-helix) generation for all crosslink types required first, and then the mixing script with force field specified but topology generation off.

i.e. sample workflow in colbuilder/colbuilder-dev environment
1. triple_helix_HLKNL.yaml (sequence generate w/ crosslink HLKNL at desired position, no geometry or topology generator, debug: true)
2. triple_helix_PYD.yaml (repeat for crosslink PYD, debug: true)
3. repeat _ad infinitum_ for however many crosslink types you want
3. mix_geometry.yaml 
 - mix_bool with ratio of two/more input pdb from above
 - set contact_distance and fibril_length as desired
 - set force_field: "martini3"
 - sequence/geometry/topology_generator: false
4. python topology_from_mix.py mix_geometry.yaml
 - this *should* pull all requisite files together and build geometry
