include $(EGS_CONFIG)
include $(SPEC_DIR)egspp1.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

# Specify the name of the user code.
# The name of the executable is determined from this variable.
#
USER_CODE = egs_mird

# The following can be used to add user macros and mortran subroutines.
# The file(s) specified here are added after egsnrc.macros, machine.macros
# and egs_c_interface2.macros but before any files that have
# executable code.
#
EGSPP_USER_MACROS = egs_mird.macros

# Specify from which base class this application is being derived.
# This has the effect of automatically compiling the base application
# class and including it into the list of object files.
#
EGS_BASE_APPLICATION = egs_advanced_application

# Specify the set of mortran sources used.
# Here we simply use the standard set defined in $HEN_HOUSE/specs/egspp1.spec
#
CPP_SOURCES = $(C_ADVANCED_SOURCES)

# Specify here other header files that your user code depends upon.
#
other_dep_user_code = $(ABS_EGSPP)egs_scoring.h $(ABS_EGSPP)geometry$(DSEP)egs_glib$(DSEP)egs_glib.h $(ABS_EGSPP)geometry$(DSEP)egs_nd_geometry$(DSEP)egs_nd_geometry.h

include $(HEN_HOUSE)makefiles$(DSEP)cpp_makefile
