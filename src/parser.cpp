#include "parser.h"
#include "utils.h"

#include <string.h>
#include <stdio.h>
#include <limits.h>

static int searcharg(const char* arg, const char* elem)
{
	return strncasecmp(arg, elem, strlen(elem));
}

/*void readint(const char* s, int* r)
{	unsigned l = UINT_MAX;
	if (sscanf(s, "%i%n", r, &l) != 1 || l != strlen(s))
		die(42, "Integer value expected, found %s", s);
}*/
void OptionDesc::ParseArg(const char* value, unsigned* param) const
{	bool ex = *value == '^';
	if (ex)
		++value;
	unsigned l = UINT_MAX;
	if (sscanf(value, "%u%n", param, &l) != 1 || l != strlen(value))
		die(42, "Unsigned integer value expected for option %s, found '%s'.", Option, value);
	if (ex)
	{	if (*param >= sizeof(unsigned)*CHAR_BIT)
			die(42, "Power of 2 constant ^%u exceeds the domain of unsigned integer (option %s).", *param, Option);
		*param = 1 << *param;
	}
}
void OptionDesc::ParseArg(const char* value, unsigned* param, unsigned min, unsigned max) const
{	ParseArg(value, param);
	if (*param < min || *param > max)
		die(42, "Value %u of option %s is out of range [%u,%u].", *param, Option, min, max);
}

void OptionDesc::ParseArg(const char* s, double* r) const
{	unsigned l = UINT_MAX;
	if (sscanf(s, "%lf%n", r, &l) != 1 || l != strlen(s))
		die(42, "Floating point value expected for option %s, found '%s'.", Option, s);
}

void OptionDesc::ParseArg(const char* s, bool* r) const
{	if (*s == 0 || strcasecmp(s, "toggle") == 0)
		*r = !*r;
	else if (strcmp(s, "1") == 0 || strcasecmp(s, "true") == 0 || strcasecmp(s, "yes") == 0 || strcasecmp(s, "on") == 0)
		*r = true;
	else if (strcmp(s, "0") == 0 || strcasecmp(s, "false") == 0 || strcasecmp(s, "no") == 0 || strcasecmp(s, "off") == 0)
		*r = false;
	else
		die(42, "Invalid boolean value '%s' for option %s.", s, Option);
}

void OptionDesc::ParseArg(const char* value, const char** param) const
{	*param = strdup(value); // TODO: memory leak
}

void Parser::HandleArg(const char* arg)
{
	if (arg[0] == '@' || arg[0] == '<')
	{	// indirect file
		FILEguard cf(arg + 1, "r");
		while (!feof(cf))
		{	char buffer[1024];
			fgets(buffer, sizeof buffer, cf);
			size_t l = strlen(buffer);
			if (l >= sizeof buffer - 1 && buffer[sizeof buffer - 2] != '\n')
				die(39, "Line in command file exceeds 1023 characters.");
			// strip whitespaces
			while (strchr(" \t\r\n", buffer[--l]) != NULL)
				buffer[l] = 0;
			if (buffer[0] == 0)
				continue; // skip empty lines
			const char* ap = buffer;
			while (strchr(" \t\r\n", *ap) != NULL)
				++ap;
			if (ap[0] == '#')
				continue; // skip comments
			HandleArg(ap);
		}
		return;
	}
	OptionDesc* ap = (OptionDesc*)bsearch(arg, ArgMap, ArgMapSize, sizeof *ArgMap, (int (*)(const void*, const void*))&searcharg);
	if (ap == NULL)
		die(44, "Illegal option %s.", arg);
	arg += strlen(ap->Option);
	if (*arg == '=')
		++arg;
	(ap->Handler)(*ap, arg);
}

void Parser::PrintHelp() const
{

}
