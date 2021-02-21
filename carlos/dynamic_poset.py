import poset
import dynamic_class
import poset_storage
from importlib import reload
reload(poset)
reload(dynamic_class)
reload(poset_storage)

Poset = dynamic_class.dynamic_class(poset.Poset)
PosetStorage = dynamic_class.dynamic_class(poset_storage.PosetStorage)
