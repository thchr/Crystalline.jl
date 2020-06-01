import bpy
import os
import sys

argv = sys.argv
argv = argv[argv.index("--") + 1:]  # get all args after "--"
# we assume by convention the following input ordering: 
#   argv[0]: savename (a string)
#   argv[1]: render engine (a string: 'eevee' or 'cycles')
#   argv[2]: TODO: Resolution?

dirpath=os.path.dirname(os.path.abspath(__file__))
#dirpath=os.path.dirname(os.path.abspath(bpy.data.filepath))

# load existing "template" .blend file, with the materials/lighting we need
bpy.ops.wm.open_mainfile(filepath=os.path.join(dirpath, "bpy-template.blend"))

# set render engine
render_engine = 'BLENDER_EEVEE' if len(argv) == 1 else argv[1].upper()
if render_engine == 'EEVEE': render_engine = 'BLENDER_'+render_engine
bpy.data.scenes[0].render.engine = render_engine

# import the obj files
bpy.ops.import_scene.obj(filepath=os.path.join(dirpath, "tmp", "isocaps.obj"))
bpy.ops.import_scene.obj(filepath=os.path.join(dirpath, "tmp", "isosurface.obj"))
bpy.ops.import_scene.obj(filepath=os.path.join(dirpath, "tmp", "unitcell.obj"))
# bpy.data.meshes["isosurface"].use_auto_smooth

# set their material properties
bpy.data.objects['isocaps'].material_slots[0].material = bpy.data.materials['Caps']
bpy.data.objects['isosurface'].material_slots[0].material = bpy.data.materials['Isosurface']

# do some tricks to get the unit cell to render as a (thin!) curve
bpy.data.objects['unitcell'].select_set(True) # select object
bpy.context.view_layer.objects.active = bpy.data.objects['unitcell'] # select^2: so fucked, see https://blender.stackexchange.com/questions/38618/selecting-an-object-via-scripting-in-2020
bpy.ops.object.editmode_toggle() # bug in .obj with lines; so fucked, see https://blender.stackexchange.com/q/153149/78703
bpy.ops.object.editmode_toggle()
bpy.ops.object.convert(target='CURVE') # convert lines to curve
bpy.data.curves['unitcell'].bevel_object = bpy.data.objects['Circle'] # "solidify" using auxilary circle curve
bpy.ops.object.material_slot_add() # add a material slot to curve
bpy.data.curves['unitcell'].materials[0] = bpy.data.materials['Unitcell'] # set to transparent material
bpy.data.collections['unitcell-freestyle'].objects.link(bpy.context.active_object) # add active object to collection 'uc'
bpy.ops.transform.resize(value=(0.999, 0.999, 0.999)) # shrink it a little, so we don't get overlapping freestyle lines

# render and save (for some reason, have to specify absolute path)
bpy.context.scene.render.filepath = os.path.join(dirpath, argv[0])
bpy.ops.render.render(write_still = True)
