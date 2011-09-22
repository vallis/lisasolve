% Whitepaper on good practice for correctness and reliability in gravitational-wave data analysis
% Michele Vallisneri
% 2011-09-21

# Coding it

## General remarks, veering on philosophy

### Writing code to communicate, not implement

We (scientists) are used to thinking of computer code as providing the *implementation* of numerical or algorithmic idea; therefore, we measure the success of a coding endeavor by the correctness, ease of use, and speed (hopefully, in this order) of the resulting computer program. I would argue instead that the practice of scientific programming would be helped by a shift of perspective to viewing computer code as embodying the successful *communication* of scientific recipes. The recipients of such communication include our collaborators (including our future selves) who will inspect and re-use the code; other users in the scientific community, if (better, when!) the code is made publicly available; referees at scientific journals, if (when!) an article is accompanied by the code used to produce its numbers and figures. From this viewpoint, using the code to actually run a calculation is a happy byproduct of communicating scientific ideas.

Such a principle immediately brings to mind favorite tools for embedding comments and documentation inside code, such as the javadoc-inspired [doxygen] for several languages, docstrings and [sphinx] for Python, and the more radical (but somewhat [flawed]) path of [literate programming]. Mathematica's approach is closer to the reverse – embedding code inside descriptive text – but the resulting notebooks are not always directly usable as modular code, or in a batch mode.

However, when I talk of *communicative code*, I am thinking about the code itself rather than comments and documentation (which are still very important, no question there). Thus, I am partial to computer languages with *expressive* syntax that can be read immediately as the mathematics or algorithms that they embody. Indeed, I shall be more audacious, and suggest that *beauty* can serve as the first line of defense from complexity and from bugs. Beautiful code manifests the parallelism between a program and its mathematical and physical content; ugly code hides it and provide a breeding ground for inconsistencies and inefficiencies.

True, there is certainly an element of taste and habit in what constitutes beautiful code. Andy Oram and Greg Wilson's [anthology], indeed entitled _Beautiful Code_ [@OramWilson2007], collects many examples in different computer languages, some of which look positively frightening to my eye. For instance, a chapter on linear-algebra routines in LAPACK showcases rather forbidding Fortran code with little trace of the elegance and synthesis of matrix notation (MATLAB does better in that regard). However, the beauty of code should be judged in the context of its purpose and function, and, as it were, through the eyes of the code's users and perusers of the code. I am reminded of what [Gerald Jay Sussman] replied to my comment, certainly a FAQ, that his Lisp code was overrun with parentheses: "Ah, but I don't even see those anymore."

I should then weaken my claim somewhat, and posit that in our relevant context (working collaboratively within a moderately sized group to write software that analyzes gravitational-wave data), it *is* possible to find a combination of computer language and coding style that, although perhaps difficult to define, *would* look beautiful and transparent to most collaborators.[^supreme] There are certainly examples of such beautiful code in the current LAL codebase. Instead of trying to propose a definition, though, let me cite the first nine mantras in Tim Peters' *Zen of Python*[^zen]:

1. Beautiful is better than ugly.
2. Explicit is better than implicit.
3. Simple is better than complex.
4. Complex is better than complicated.
5. Flat is better than nested.
6. Sparse is better than dense.
7. Readability counts.
8. Special cases aren't special enough to break the rules.
9. Although practicality beats purity.

Mantras 1, 3, and 7 seem especially relevant to this discussion: beauty, simplicity, and readability allow the clear and economical expression of mathematical concepts; they provide a transparent representation of the flow of information; and they help build a first line of defense against bugs (at least unsubtle ones). I will discuss [below](#scripting-languages) how Python, which was built according to these principles, and which has a user community that actively encourages their adoption, can *feel like home*[^waits] for a scientist concerned about correctness, usability, and efficiency. Nevertheless, the same principles can be applied to other languages that may very well be more appropriate to some applications.

#### Resources and diversions

* David Gelernter, more eloquently than I can, highlights the [importance of beauty](gelernter) in computer languages and programs in his 1998 book [@Gelernter1998].

* In his essay on "Treating Code As An Essay" [@OramWilson2007], Ruby's creator Yukihiro Matsumoto lists the qualities of beautiful code as brevity, familiarity, simplicity, flexibility, and balance. (Ruby is the Japanese response to Python, and it takes the zen to a maximum. It is also the language of choice in Piet Hut and Jun Makino's [Art of Computational Science](acs) magnum opus.)

* There is of course much to say about the role of beauty in *physics* (and especially in theoretical physics), where it has been a guiding principle for the development of more and more general, economical, and powerful theories. Chandrasekhar may have said it best in his [essay](chandra) and book [@Chandrasekhar1987].

* For beauty in the *graphical display of information*, see *Plotting it* in this very report.

[doxygen]: http://www.doxygen.org "Doxygen documentation system"
[sphinx]: http://sphinx.pocoo.org "Sphinx Python documentation generator"
[flawed]: http://software-carpentry.org/2011/03/4069 "Greg Wilson's opinion piece on Literate Programming"
[literate programming]: http://en.wikipedia.org/wiki/Literate_programming "Literate Programming in wikipedia"
[anthology]: http://shop.oreilly.com/product/9780596510046.do "Oram and Wilson: Beautiful Code"
[lapack]: http://netlib.org/utk/people/JackDongarra/PAPERS/beautiful-code.pdf "Jack Dongarra on LAPACK's beauty"
[Gerald Jay Sussman]: http://groups.csail.mit.edu/mac/users/gjs "Gerald Jay Sussman's home page"
[chandra]: http://history.fnal.gov/GoldenBooks/gb_chandrasekhar.html "Subrahmanyan Chandrasekhar: Beauty and the Quest for Beauty in Science"
[gelernter]: http://www.theatlantic.com/past/docs/unbound/digicult/dg4.htm "David Gelernter: excerpt from Machine Beauty: Elegance and the Heart of Technology"
[acs]: http://www.artcompsci.org "Hut and Makino's Art of Computational Science"

[^supreme]: Sort of a negative Miller test, isn't it?
[^zen]: `>>> import this`
[^waits]: Although "forgive me pretty baby but I always take the long way home" (Norah Jones/Tom Waits).

### Computer code as a non-native language for scientists

While my exhortation that we seek beautiful code has an aspirational character, in this section I must strike a more cautionary note: we should always remember that, for must of us, *computer code is a non-native language* that we will always write, as it were, with an accent, no matter how much we enjoy doing so. When I asked reliable-software expert Gerard Holzmann what he thought of scientists writing software, he replied "What would you think of programmers writing equations for you?" I certainly wouldn't go as far as that: anecdotal evidence suggests[^wrong] that embedding programmers in research groups, even when financially feasible, is not too efficient because, in a sense, scientific software is closer to the equations than to the machine.[^sociology]

There is a more fundamental objection to Holzmann's reply, which is that  software-development theory, which presumably informs the M.O. of professional programmers, [does not really apply](scientific programming) to most of scientific programming. The standard cycle of software development addressed by formal methodologies (say, [Agile] or [Extreme] Programming) could be caricatured as formulating specifications with a *customer* and then working through milestones and releases that satisfy them. But scientists are often their own customers, seldom write specifications (indeed, they would hardly know what to write before they have started using the software), rarely bother with releases, but rather continually exercise their evolving code until its science output makes sense to them.

In these conditions, it may not make sense to prescribe a methodology for scientific programming, but rather we should encourage scientists to adopt *best practices* that can help them move through the non-native world of computation. Such best practices would include using (in no particular order) modern languages, version control, operating-system facilities, testing, object-oriented techniques, databases, and so on. This is the spirit of Greg Wilson's excellent [Software Carpentry](http://software-carpentry.org) course, which promises to "make scientists and engineers dramatically more productive by teaching them fundamental computational skills." I would support this appeal to self interest with an exhortation toward common good: the correctness, reliability, and reproducibility of computational research endeavors is increasingly influential and important for society at large, which devotes precious resources to funding researchers. Thus, we *need* scientists to embrace good computational practices that can make them more productive *and* make their results more reliable.[^sociology2]

[holzmann]: http://lars-lab.jpl.nasa.gov "JPL Laboratory for Reliable Software"
[scientific programming]: http://software-carpentry.org/articles/how-scientists-use-computers-2009.pdf "Hannay et al.: How Do Scientists Develop and Use Scientific Software?"
[agile]: http://en.wikipedia.org/wiki/Agile_software_development
[extreme]: http://en.wikipedia.org/wiki/Extreme_programming

[^wrong]: So I may be *completely* wrong, so I would like to see data on this!
[^sociology]: Not that turning graduate students and postdocs into programmers is that great, either. In this age, computational... drudgery should be part of the balanced diet of every practicing scientist, but it should not, for instance, dominate the day-to-day work of young physicists, not least because it is not a valued contribution career-wise. The reader should gratefully notice that I confining my sociological gripes to a footnote, at least for the moment.
[^sociology2]: Sociological footnote no. 2. In practice, a few things need to happen to make such exhortations effective: undergraduate and postgraduate curricula need to make space for computation and software carpentry; journals and funding agencies must require researchers to apply the same care and high standards to software that they devote to article submissions; hiring and tenure committees must attribute value to scientific-software contributions.

#### Resources and diversions

* Greg Wilson's [Software Carpentry](http://software-carpentry.org).
* A somewhat orthogonal topic of teaching to software carpentry is [Computational Thinking](http://www.cs.cmu.edu/~CompThink). Also [Jon Udell's take](http://www.slideshare.net/judell/computational-thinking), and [Natural Programming](http://www.cs.cmu.edu/~NatProg/index.html) projects.
* Google's [Code University](http://code.google.com/edu).
* "Computational-X" was the third paradigm... "X-informatics" is the [fourth](http://research.microsoft.com/en-us/collaboration/fourthparadigm/). E.g., [AstroInformatics](http://www.astroinformatics2010.org).

### An earnest pleading for open code (and open data)

#### Resources and diversions

* Reproducible research. *Computing in Science and Engineering* 11(1), 2009.
* Victoria Stodden: [Trust Your Science? Open Your Data and Code](http://magazine.amstat.org/blog/2011/07/01/trust-your-science).
* [Panton Principles](http://pantonprinciples.org) for open data in science.

### More

#### Resources and diversions

## A very simple taxonomy of languages and their use in GW data analysis

### Compiled languages

#### Resources

* LAL and XLAL: LIGO Scientific Collaboration Algorithm Library [Specification and Style Guide](https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=T990030).
* JPL's Institutional [Coding Standard](http://lars-lab.jpl.nasa.gov/JPL_Coding_Standard_C.pdf) for the C Programming Language.

### Scripting languages

Python can help because:

* It is expressive (Python/C = 6) and introspective.
* It is object-oriented, but casually so.
* It provides high-level mathematical–numerical objects and libraries.
* It can transparently wrap high-efficiency compiled code.

#### Resources

* [Python: Batteries Included](http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=4160244). *Computing in Science and Engineering* 9(3), 2007.

### Integrated numerical and symbolic environments

## References
