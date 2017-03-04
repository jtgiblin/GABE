GABE
====

Grid And Bubble Evolver - Emscripten branch

GABE gutted down to the basics (no fft's, no file I/O) for use in a browser.
To compile, you will need to install [https://kripken.github.io/emscripten-site/](emscripten).

For now, just compile using something like
```
em++ -std=c++11 g2functions.cpp  g2init.cpp  g2main.cpp  g2model.cpp  g2output.cpp
```
