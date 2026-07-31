# Long chemistry string regression fixture

This mechanism exercises the two legacy truncation boundaries that used to
break otherwise valid CHEMKIN input:

- six species identifiers occupy all 18 columns reserved by NASA-7;
- the `SPECIES` declaration and reaction record are both longer than 80
  characters.

All six long-name species reuse the same H2 elemental composition and
thermodynamic polynomial.  The single three-to-three reaction is therefore
element-balanced; its purpose is parser and cklink round-trip coverage, not
physical kinetics.
