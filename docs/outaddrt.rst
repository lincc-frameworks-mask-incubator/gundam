.. _outaddrt:

=======================================
Additional Runtime Information in Ouput
=======================================

Some useful information is added to the **par** dictionary returned within the
**counts** dictionary. These are mostly to keep a record of runtime parameters,
such as if the run was parallel, if it wrote ouput files, check the boundaries 
of the survey, etc.
    
+----------------+--------------------------------------------------------+
| Key(s)         | Description                                            |
+================+========================================================+
| para           | Wheter to run in parallel or serial, i.e single core   |
+----------------+--------------------------------------------------------+
| write          | Write files was requested or not                       |
+----------------+--------------------------------------------------------+
| plot           | Plot was requested or not                              |
+----------------+--------------------------------------------------------+
| sbound         | Survey boundaries as a list (ramin,ramax,decmin,       |
|                | decmax,dcmin,dcmax), as returned for example by        |
|                | :func:`gundam.bound3d`                                 |
+----------------+--------------------------------------------------------+
