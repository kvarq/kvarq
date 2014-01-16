
import os.path

from kvarq.genes import Genome

ancestor_path = os.path.join(os.path.dirname(__file__), os.path.pardir, 'MTB_ancestor_reference.bases')
ancestor = Genome(ancestor_path, 'MTB ancestor')
# win32 GIT checkout can add '\r'
assert ancestor.size == 4411533 or ancestor.size == 4411534

