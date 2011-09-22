all: html/coding.html html/running.html html/sharing.html html/plotting.html

html/%.html: src/%.md src/refs.mods
	pandoc --bibliography src/refs.mods --standalone -f markdown -t html $< | tidy -utf8 -i | cat > $@

src/refs.mods: src/refs.bib
	bib2xml $< > $@
