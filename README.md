# edf2cfs
A high-performance multi-threaded C++ CLI application to convert EDF (European Data Format) files to Compressed Feature Set (CFS) format. The CFS format is used by the Z3Score sleep scoring system (https://z3score.com). Instead of using polysomnography data in European Data Format (EDF, https://en.wikipedia.org/wiki/European_Data_Format), the Z3Score system uses CFS files. CFS files are on an average 17X smaller than corresponding EDF files. This reduces data overhead significantly. The format does not allow any user identifiable information ensuring anonymity.

Patents pending (c)-2017 Amiya Patanaik amiyain@gmail.com

This application has miltiple dependencies including 

sigpack: C++ signal processing library http://sigpack.sourceforge.net/ 

tclap: Templatized C++ Command Line Parser Library, 

BOOST: C++ libraries 

armadillo: C++ linear algebra library

zlib:  DEFLATE compression algorithm

If you have all the dependencies installed you can issue the command

```sh
	make
```

edf2cfs is the compiled executable for Mac-OS, to find out the usage

```
	./edf2cfs -h
```

```
USAGE: 

   ./edf2cfs  [-l] [-o] [-q] [-d <EDF Directory>] [-z <ER-A1 Channel
              Label>] [-x <EL-A2 Channel Label>] [-b <C4-A1 Channel Label>]
              [-a <C3-A2 Channel Label>] [--] [--version] [-h] <List of EDF
              files> ...


Where: 

   -l,  --log
     save log

   -o,  --overwrite
     over write files

   -q,  --quiet
     silent mode

   -d <EDF Directory>,  --dir <EDF Directory>
     EDF Directory

   -z <ER-A1 Channel Label>,  --er <ER-A1 Channel Label>
     ER-A1 Channel Label

   -x <EL-A2 Channel Label>,  --el <EL-A2 Channel Label>
     EL-A2 Channel Label

   -b <C4-A1 Channel Label>,  --c4 <C4-A1 Channel Label>
     C4-A1 Channel Label

   -a <C3-A2 Channel Label>,  --c3 <C3-A2 Channel Label>
     C3-A2 Channel Label

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <List of EDF files>  (accepted multiple times)
     List of filename


   Usage: ./edf2cfs -a C3A2 -b C4A1 -x ELA2 -z ERA1 -q -o -l -d edfDir
   filename1.edf filename2.edf ... filenameN.edf

   If no channels are given, then a selection menu will be shown.

   Use -d to provide a directory path with EDF files, -q to supress output,
   -o to overwrite and -l to save log.

```

Both batch and single file conversion is supported

License
----

GPL V3
