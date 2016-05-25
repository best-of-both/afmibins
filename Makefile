#########################################################################
#                                                                       #
# Immersed Boundary Incompressible Navier-Stokes solver                 #
#                                                                       #
# Copyright (C) 2016  Andrew Kassen <atkassen@gmail.com>                #
#                                                                       #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#                                                                       #
#########################################################################

rootdir = $(shell git rev-parse --show-cdup)
need_modules = 
include $(rootdir)Makefile.in

OBJDIR = obj
LIBDIR = lib
DEPDIR = deps
BINDIR = bin
INCDIR = include

.PHONY: clean distclean

clean:
	-@rm -r $(OBJDIR) $(LIBDIR) 2> /dev/null || true

distclean:
	-@rm -r $(OBJDIR) .$(INCDIR) $(DEPDIR) $(LIBDIR) $(BINDIR) 2> /dev/null || true
