This repo adheres to a fairly strict coding style.  Though we aim to be rigorous in following this style, if you find it impeding your productivity, we encourage you to bring forth evidence as to why.  The important thing is that we all adhere to the same standards so we may operate as one!  Further, our style is borrowed largely, if not entirely, from the coding community at large.  As such, all of these practices should bring smiles to the faces of most collaborators you will encounter in future endeavors.

Code
-----

Basically, this is (PEP 8)[http://www.python.org/dev/peps/pep-0008/]
but we also make some specific choices that are given as options in pep 8.

- Classes are all upper case and camel case: ClassName
- methods are all lowercase with underscores: method_name() or _private_method_name()

pylint - We also look to follow (pylint)[http://www.pylint.org/] recomendations.
    This is nice because you can install the sublimelinter package and it will check your code for you.  This makes adhering to pep 8 and lint rules really easy because it highlights any problems.

line width:
    - 80 chars is great
    - 100 chars is fine
    - 120 chars is prob too long.  figure a way to break it up and get it close to 100

variable names - these should also be all lower case with underscores (variable_name) except for constants.
    Additionally, abbreviations should be mostly avoided except for very short lived variables and extremely common abbreviations specific to this project.  When it doubt, spell it out. This can make lines kind of long, but that's one reason we are lax about line length.  Also, if we ultimately want to shorten stuff, it will be very easy to refactor explicit variable names.
