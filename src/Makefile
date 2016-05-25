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
currdir = $(shell pwd)
root = $(realpath $(currdir)/$(rootdir))
reldir = $(strip $(subst $(root)/,,$(currdir)))

OBJDIR ?= obj
DEPDIR ?= deps
INCDIR ?= include
BINDIR ?= bin
SRCDIR ?= src
mangled = $(subst //,/,$(subst $(SRCDIR),,$(reldir)))
ifeq "$(mangled)" ""
mangled = root
endif

need_modules = 
makefile ?= Makefile

submodules = $(shell ls -F */$(makefile) 2> /dev/null)

green = "\033[32;1m"
red = "\033[31;1m"
yellow = "\033[33;1m"
normal = "\033[0;0m"
depth = [$(yellow)$(mangled)$(normal)]

.PHONY: all test clean distclean FORCE

all: $(submodules)

test: $(submodules)

clean: $(submodules)

distclean: $(submodules)

%/Makefile: FORCE
	@echo $(depth) $(green)make$(normal) $* $(MAKECMDGOALS)
	@$(MAKE) --silent --directory=$* $(MAKECMDGOALS)

FORCE:

