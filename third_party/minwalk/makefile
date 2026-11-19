CXX := g++
CXXFLAGS := -Wall -O3 -std=c++11
ifeq ($(OS),Windows_NT)
    CXXFLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
        CXXFLAGS += -D AMD64
    endif
    ifeq ($(PROCESSOR_ARCHITECTURE),x86)
        CXXFLAGS += -D IA32
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CXXFLAGS += -D LINUX
        LIBS = -lpthread -lglut -lGL -lGLU
    endif
    ifeq ($(UNAME_S),Darwin)
        CXXFLAGS += -D OSX
        LIBS = -framework GLUT -framework OpenGL
    endif
endif

CXXFLAGSGUI := $(CXXFLAGS) -Wno-deprecated -D GRAPH

SRCDIR := src
BUILDDIR := build
DISTDIR := dist

TARGET := hpfold
TARGETGUI := hpfoldgui

OBJECTS := $(addprefix $(BUILDDIR)/,Grid.o HPFolding.o main.o)
OBJECTSTEST := $(addprefix $(BUILDDIR)/,Grid.o HPFolding.o test.o)

OBJECTSGUI := $(addprefix $(BUILDDIR)/,Grid.o HPFolding.o gHPGUI.o gmain.o)



$(TARGET): $(OBJECTS)
	@echo " Linking..."
	mkdir -p $(DISTDIR)
	@echo " $(CXX) $^ -o $(DISTDIR)/$(TARGET) "; $(CXX) $^ -o $(DISTDIR)/$(TARGET)

$(TARGETGUI): $(OBJECTSGUI)
	@echo " Linking..."
	mkdir -p $(DISTDIR)
	@echo " $(CXX) $^ -o $(DISTDIR)/$(TARGETGUI) $(LIBS) "; $(CXX) $^ -o $(DISTDIR)/$(TARGETGUI) $(LIBS)

tests: $(OBJECTSTEST)
	@echo " Linking..."
	mkdir -p $(DISTDIR)
	@echo " $(CXX) $^ -o $(DISTDIR)/hptest "; $(CXX) $^ -o $(DISTDIR)/hptest


$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CXXFLAGS) -c -o $@ $<"; $(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILDDIR)/%.o: tests/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CXXFLAGS) -c -o $@ $<"; $(CXX) $(CXXFLAGS) -c -o $@ $<


$(BUILDDIR)/g%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CXXFLAGSGUI) -c -o $@ $<"; $(CXX) $(CXXFLAGSGUI) -c -o $@ $<




clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(DISTDIR)"; $(RM) -r $(BUILDDIR) $(DISTDIR)
