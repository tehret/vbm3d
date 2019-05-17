#ifndef __CMD_OPTION_H__
#define __CMD_OPTION_H__

#include <cstdarg>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

//
// This is a copy paste of CImg option functions, so that they can be used
// independently of CImg. To avoid name colisions, cimg_ prefixes have
// been replaced by clo_.
//



// Define the program usage, and retrieve command line arguments.
//
#define clo_usage(usage) \
  cmd_opt::option((char*)0,argc,argv,(char*)0,usage)
#define clo_help(str) \
  cmd_opt::option((char*)0,argc,argv,str,(char*)0)
#define clo_option(name,defaut,usage) \
  cmd_opt::option(name,argc,argv,defaut,usage)
#define clo_argument(pos) \
  cmd_opt::argument(pos,argc,argv)
#define clo_argument1(pos,s0) \
  cmd_opt::argument(pos,argc,argv,1,s0)
#define clo_argument2(pos,s0,s1) \
  cmd_opt::argument(pos,argc,argv,2,s0,s1)
#define clo_argument3(pos,s0,s1,s2) \
  cmd_opt::argument(pos,argc,argv,3,s0,s1,s2)
#define clo_argument4(pos,s0,s1,s2,s3) \
  cmd_opt::argument(pos,argc,argv,4,s0,s1,s2,s3)
#define clo_argument5(pos,s0,s1,s2,s3,s4) \
  cmd_opt::argument(pos,argc,argv,5,s0,s1,s2,s3,s4)
#define clo_argument6(pos,s0,s1,s2,s3,s4,s5) \
  cmd_opt::argument(pos,argc,argv,6,s0,s1,s2,s3,s4,s5)
#define clo_argument7(pos,s0,s1,s2,s3,s4,s5,s6) \
  cmd_opt::argument(pos,argc,argv,7,s0,s1,s2,s3,s4,s5,s6)
#define clo_argument8(pos,s0,s1,s2,s3,s4,s5,s6,s7) \
  cmd_opt::argument(pos,argc,argv,8,s0,s1,s2,s3,s4,s5,s6,s7)
#define clo_argument9(pos,s0,s1,s2,s3,s4,s5,s6,s7,s8) \
  cmd_opt::argument(pos,argc,argv,9,s0,s1,s2,s3,s4,s5,s6,s7,s8)

namespace cmd_opt {

	const char t_normal[] = { 0x1b,'[','0',';','0',';','0','m','\0' };
	const char t_red[]    = { 0x1b,'[','4',';','3','1',';','5','9','m','\0' };
	const char t_bold[]   = { 0x1b,'[','1','m','\0' };
	const char t_purple[] = { 0x1b,'[','0',';','3','5',';','5','9','m','\0' };
	const char t_green[]  = { 0x1b,'[','0',';','3','2',';','5','9','m','\0' };

	//! Remove the 'case' of an ASCII character.
	inline char uncase(const char x) {
		return (char)((x<'A'||x>'Z')?x:x-'A'+'a');
	}

	//! Remove the 'case' of a C string.
	/**
	  Acts in-place.
	 **/
	inline void uncase(char *const string) {
		if (string) for (char *ptr = string; *ptr; ++ptr) *ptr = uncase(*ptr);
	}

	//! Read a float number from a C-string.
	/**
	  \note This function is quite similar to <tt>std::atof()</tt>,
	  but that it allows the retrieval of fractions as in "1/2".
	 **/
	inline float atof(const char *str) {
		float x = 0,y = 1;
		if (!str) return 0; else { sscanf(str,"%g/%g",&x,&y); return x/y; }
	}

	//! Compute the length of a C-string.
	/**
	  \note This function is similar to <tt>std::strlen()</tt>
	  and is here because some old compilers do not
	  define the <tt>std::</tt> version.
	 **/
	inline int strlen(const char *s) {
		if (s) { int k; for (k = 0; s[k]; ++k) ; return k; }
		return -1;
	}

	//! Compare the first \p n characters of two C-strings.
	/**
	  \note This function is similar to <tt>std::strncmp()</tt>
	  and is here because some old compilers do not
	  define the <tt>std::</tt> version.
	 **/
	inline int strncmp(const char *s1, const char *s2, const int l) {
		if (s1 && s2) { int n = 0; for (int k = 0; k<l; ++k) n+=abs(s1[k] - s2[k]); return n; }
		return 0;
	}

	//! Compare the first \p n characters of two C-strings, ignoring the case.
	/**
	  \note This function is similar to <tt>std::strncasecmp()</tt>
	  and is here because some old compilers do not
	  define the <tt>std::</tt> version.
	 **/
	inline int strncasecmp(const char *s1, const char *s2, const int l) {
		if (s1 && s2) { int n = 0; for (int k = 0; k<l; ++k) n+=abs(uncase(s1[k])-uncase(s2[k])); return n; }
		return 0;
	}

	//! Compare two C-strings.
	/**
	  \note This function is similar to <tt>std::strcmp()</tt>
	  and is here because some old compilers do not
	  define the <tt>std::</tt> version.
	 **/
	inline int strcmp(const char *s1, const char *s2) {
		const int l1 = strlen(s1), l2 = strlen(s2);
		return strncmp(s1,s2,1+(l1<l2?l1:l2));
	}

	//! Compare two C-strings, ignoring the case.
	/**
	  \note This function is similar to <tt>std::strcasecmp()</tt>
	  and is here because some old compilers do not
	  define the <tt>std::</tt> version.
	 **/
	inline int strcasecmp(const char *s1, const char *s2) {
		const int l1 = strlen(s1), l2 = strlen(s2);
		return strncasecmp(s1,s2,1+(l1<l2?l1:l2));
	}

	//! Find a character in a C-string.
	inline int strfind(const char *s, const char c) {
		if (s) {
			int l; for (l=strlen(s); l>=0 && s[l]!=c; --l) ;
			return l;
		}
		return -1;
	}

	//! Compute the basename of a filename.
	inline const char* basename(const char *s)  {
		return ("this is not windows")?(s?s+1+strfind(s,'/'):0):(s?s+1+strfind(s,'\\'):0);
	}



	inline const char* option(const char *const name, const int argc, char **argv,
			const char *defaut, const char *const usage=0) {
		static bool first = true, visu = false;
		const char *res = 0;
		if (first) {
			first=false;
			visu = (option("-h",argc,argv,(char*)0)!=0);
			visu |= (option("-help",argc,argv,(char*)0)!=0);
			visu |= (option("--help",argc,argv,(char*)0)!=0);
		}
		if (!name && visu) {
			if (usage) {
				fprintf(stderr,"\n %s%s%s", t_red,basename(argv[0]),t_normal);
				fprintf(stderr," : %s",usage);
				fprintf(stderr," (%s, %s)\n\n",__DATE__,__TIME__);
			}
			if (defaut) fprintf(stderr,"%s\n",defaut);
		}
		if (name) {
			if (argc>0) {
				int k = 0;
				while (k<argc && strcmp(argv[k],name)) ++k;
				res = (k++==argc?defaut:(k==argc?argv[--k]:argv[k]));
			} else res = defaut;
			if (visu && usage) fprintf(stderr,"    %s%-8s%s = %-12s : %s%s%s\n",
					t_bold,name,t_normal,res?res:"0",t_purple,usage,t_normal);
		}
		return res;
	}

	inline bool option(const char *const name, const int argc, char **argv,
			const bool defaut, const char *const usage=0) {
		const char *s = option(name,argc,argv,(char*)0);
		const bool res = s?(strcasecmp(s,"false") && strcasecmp(s,"off") && strcasecmp(s,"0")):defaut;
		option(name,0,0,res?"true":"false",usage);
		return res;
	}

	inline int option(const char *const name, const int argc, char **argv,
			const int defaut, const char *const usage=0) {
		const char *s = option(name,argc,argv,(char*)0);
		const int res = s?atoi(s):defaut;
		char tmp[256];
		sprintf(tmp,"%d",res);
		option(name,0,0,tmp,usage);
		return res;
	}

	inline char option(const char *const name, const int argc, char **argv,
			const char defaut, const char *const usage=0) {
		const char *s = option(name,argc,argv,(char*)0);
		const char res = s?s[0]:defaut;
		char tmp[8];
		tmp[0] = res;
		tmp[1] ='\0';
		option(name,0,0,tmp,usage);
		return res;
	}

	inline float option(const char *const name, const int argc, char **argv,
			const float defaut, const char *const usage=0) {
		const char *s = option(name,argc,argv,(char*)0);
		const float res = s?atof(s):defaut;
		char tmp[256];
		sprintf(tmp,"%g",res);
		option(name,0,0,tmp,usage);
		return res;
	}

	inline double option(const char *const name, const int argc, char **argv,
			const double defaut, const char *const usage=0) {
		const char *s = option(name,argc,argv,(char*)0);
		const double res = s?atof(s):defaut;
		char tmp[256];
		sprintf(tmp,"%g",res);
		option(name,0,0,tmp,usage);
		return res;
	}

	inline const char* argument(const unsigned int nb, const int argc, char **argv, const unsigned int nb_singles=0, ...) {
		for (int k=1, pos=0; k<argc;) {
			const char *const item = argv[k];
			bool option = (*item=='-'), single_option = false;
			if (option) {
				va_list ap;
				va_start(ap,nb_singles);
				for (unsigned int i=0; i<nb_singles; ++i) if (!strcasecmp(item,va_arg(ap,char*))) { single_option = true; break; }
				va_end(ap);
			}
			if (option) { ++k; if (!single_option) ++k; }
			else { if (pos++==(int)nb) return item; else ++k; }
		}
		return 0;
	}

} //namespace cmd_opt

#endif //__CLO_OPTION_H__

