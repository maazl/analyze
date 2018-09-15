#include "parser.h"
#include "utils.h"

#include <cstring>
#include <cstdio>
#include <climits>
#include <cmath>
#include <algorithm>

using namespace std;

static const char blanks[] = "                ";
void OptionDesc::PrintPreamble() const
{	size_t len = strnlen(Option.data(), 8);
	fwrite(Option.data(), 1, len, stderr);
	len += PrintParam();
	fputs(len < sizeof(blanks) ? blanks + len : "\t", stderr);
	fputs(Description, stderr);
}

void OptionDesc::SetBase::CheckNoVal(const char* val) const
{	if (val)
		die(42, "Option %.8s does not expect a parameter.", Option.data());
}
size_t OptionDesc::SetBase::PrintParam() const
{	return 0;
}
void OptionDesc::SetBase::DoPrint(bool active) const
{	PrintPreamble();
	if (active)
		fputs(" (selected)", stderr);
}
void OptionDesc::SetBase::DoWrite(FILE* out) const
{	fprintf(out, "%.8s\n", Option.data());
}

void OptionDesc::Opt<bool>::Parse(const char* s) const
{	if (!s || !*s)
		die(42, "Boolean option %.8s requires parameter {true|false|toggle}.", Option.data());
	if (strcasecmp(s, "toggle") == 0)
		Param = !Param;
	else if (strcmp(s, "1") == 0 || strcasecmp(s, "true") == 0 || strcasecmp(s, "yes") == 0 || strcasecmp(s, "on") == 0)
		Param = true;
	else if (strcmp(s, "0") == 0 || strcasecmp(s, "false") == 0 || strcasecmp(s, "no") == 0 || strcasecmp(s, "off") == 0)
		Param = false;
	else
		die(42, "Invalid boolean value '%s' for option %.8s.", s, Option.data());
}

void OptionDesc::Opt<unsigned>::Parse(const char* value) const
{	if (!value)
		die(42, "Option %.8s requires a value.", Option.data());
	bool ex = *value == '^';
	if (ex)
		++value;
	unsigned l = UINT_MAX;
	if (sscanf(value, "%u%n", &Param, &l) != 1 || l != strlen(value))
		die(42, "Unsigned integer value expected for option %.8s, found '%s'.", Option.data(), value);
	if (ex)
	{	if (Param >= sizeof(unsigned)*CHAR_BIT)
			die(42, "Power of 2 constant ^%u exceeds the domain of unsigned integer (option %.8s).", Param, Option.data());
		Param = 1 << Param;
	}
}

void OptionDesc::Opt<double>::Parse(const char* value) const
{	if (!value)
		die(42, "Option %.8s requires a value.", Option.data());
	int ex = *value == '^';
	if (ex)
	{	++value;
		if (*value == '/');
		{ ex = -1;
			++value;
	}	}
	unsigned l = UINT_MAX;
	if (sscanf(value, "%lf%n", &Param, &l) != 1 || l != strlen(value))
		die(42, "Floating point value expected for option %.8s, found '%s'.", Option.data(), value);
	if (ex)
		Param = pow(2, ex < 0 ? 1/Param : Param);
}

void OptionDesc::Opt<const char*>::Parse(const char* value) const
{	Param = value ? strdup(value) : NULL; // TODO: memory leak
}

void OptionDesc::LimitOpt<unsigned>::Parse(const char* val) const
{	Opt<unsigned>::Parse(val);
	if (Param < Min || Param > Max)
		die(42, "Value %u of option %.8s is out of range [%u,%u].", Param, Option.data(), Min, Max);
}
void OptionDesc::LimitOpt<double>::Parse(const char* val) const
{	Opt<double>::Parse(val);
	if (Param < Min || Param > Max)
		die(42, "Value %g of option %.8s is out of range [%g,%g].", Param, Option.data(), Min, Max);
}


void Parser::HandleArg(const char* arg)
{
	if (arg[0] == '@' || arg[0] == '<')
	{	// indirect file
		FILEguard cf(arg + 1, "r");
		while (!feof(cf))
		{	char buffer[1024];
			if (!fgets(buffer, sizeof buffer, cf))
				break;
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
	// Find argument entry
	auto op = upper_bound(ArgMap, ArgMapEnd, arg, OptionDesc::Comparer());
	if (op != ArgMap)
	{	--op;
		size_t len = strnlen(op->Ref.Option.data(), 8);
		//printf("%zu\t%s\t%s\t%c\n", len, op->Option, arg, arg[len]);
		if (strncmp(op->Ref.Option.data(), arg, len) == 0)
		{	switch (arg[len])
			{default:
				goto fail;
			 case 0: // w/o value
				arg = NULL;
				break;
			 case '=':
				arg += len + 1;
			}
			op->Ref.Parse(arg);
			return;
	}	}
 fail:
	const char* cp = strchr(arg, '=');
	int len = cp ? cp - arg : strlen(arg);
	die(44, "Illegal option %.*s.", len, arg);
}

size_t OptionDesc::Opt<bool>::PrintParam() const
{	fputs("=<bool>", stderr);
	return 7;
}
size_t OptionDesc::Opt<unsigned>::PrintParam() const
{	fputs("=<num>", stderr);
	return 6;
}
size_t OptionDesc::Opt<double>::PrintParam() const
{	fputs("=<float>", stderr);
	return 8;
}
size_t OptionDesc::Opt<const char*>::PrintParam() const
{	fputs("=<str>", stderr);
	return 8;
}

void OptionDesc::Opt<bool>::PrintValue(bool value)
{	fputs(value ? "on" : "off", stderr);
}
void OptionDesc::Opt<unsigned>::PrintValue(unsigned value)
{	fprintf(stderr, "%u", value);
}
void OptionDesc::Opt<double>::PrintValue(double value)
{	fprintf(stderr, "%g", value);
}
void OptionDesc::Opt<const char*>::PrintValue(const char* value)
{	if (value)
	{	fputc('\'', stderr);
		fputs(value, stderr);
		fputc('\'', stderr);
	} else
		fputs("<none>", stderr);
}

void OptionDesc::Opt<bool>::Print() const
{	PrintPreamble();
	fputs(" (current: ", stderr);
	PrintValue(Param);
	fputc(')', stderr);
}
void OptionDesc::Opt<unsigned>::Print() const
{	PrintPreamble();
	fputs(" (current: ", stderr);
	PrintValue(Param);
	fputc(')', stderr);
}
void OptionDesc::Opt<double>::Print() const
{	PrintPreamble();
	fputs(" (current: ", stderr);
	PrintValue(Param);
	fputc(')', stderr);
}
void OptionDesc::Opt<const char*>::Print() const
{	PrintPreamble();
	fputs(" (current: ", stderr);
	PrintValue(Param);
	fputc(')', stderr);
}

void Parser::PrintHelp() const
{	for (auto op = ArgMap; op != ArgMapEnd; ++op)
		if (op->Ref.Description)
		{	op->Ref.Print();
			fputc('\n', stderr);
		}
}

void OptionDesc::Opt<bool>::Write(FILE* out) const
{	fprintf(out, "%.8s=%c\n", Option.data(), '0' + Param);
}
void OptionDesc::Opt<unsigned>::Write(FILE* out) const
{	fprintf(out, "%.8s=%u\n", Option.data(), Param);
}
void OptionDesc::Opt<double>::Write(FILE* out) const
{	fprintf(out, "%.8s=%g\n", Option.data(), Param);
}
void OptionDesc::Opt<const char*>::Write(FILE* out) const
{	fprintf(out, Param ? "%.8s=%s\n" : "%.8s\n", Option.data(), Param);
}

void Parser::WriteConfig(FILE* out) const
{	for (auto op = ArgMap; op != ArgMapEnd; ++op)
		if (op->Ref.Description)
			op->Ref.Write(out);
}
