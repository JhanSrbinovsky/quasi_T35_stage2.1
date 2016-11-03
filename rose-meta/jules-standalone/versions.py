import rose.upgrade

from .version34_40 import *
from .version40_41 import *
from .version41_42 import *
from .version42_43 import *
from .version43_44 import *
from .version44_45 import *
from .version45_46 import *
              
class vn46_txxx(rose.upgrade.MacroUpgrade):

    """Upgrade macro from JULES by Author"""

    BEFORE_TAG = "vn4.6"
    AFTER_TAG = "vn4.6_txxx"

    def upgrade(self, config, meta_config=None):
        """Upgrade a JULES runtime app configuration."""

        # Add settings
        return config, self.reports

