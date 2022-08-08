# object files
OBJS = $(patsubst %.cpp,objs/%.o,$(notdir $(SOURCES)))

# how to make the main target
$(EXECUTABLE): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)
#	$(CC) $(LDFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: release
release: $(EXECUTABLE)

# how to compile each file
.SUFFIXES:
objs/%.o:
	$(CC) $(CFLAGS) -o $@ $<

# cleaning up
.PHONY: clean
clean:
	-rm -f objs/*.o $(MAIN_PROGRAM)

# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir objs
	-rm -f objs/depend
	$(foreach srcfile,$(SOURCES),$(CC) $(INCLUDES) -MM $(srcfile) -MT $(patsubst %.cpp,objs/%.o,$(notdir $(srcfile))) >> objs/depend;)

-include objs/depend
