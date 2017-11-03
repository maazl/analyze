#ifndef PARSER_H_
#define PARSER_H_

#include "utils.h"

#include <stdlib.h>
#include <functional>
#include <fstream>


class OptionDesc;

class Parser
{	const OptionDesc* const ArgMap;
	const size_t ArgMapSize;
 public:
	template <size_t N>
	Parser(const OptionDesc (&argmap)[N]) : ArgMap(argmap), ArgMapSize(N) {}
	void HandleArg(const char* arg);
	void PrintHelp() const;
};

struct OptionDesc
{	typedef void (*HandlerFn)(const OptionDesc& opt, const char* val);
	//typedef std::function<void(const char*)> HandlerFn;
	const char Option[8];
	const char* const Description;
	const HandlerFn Handler;
	/*template <typename ...ARGS>
	OptionDesc(const char (&option)[8], const char* description, ARGS... args)
		: Option(option), Description(description), Handler([args...](const char* val) { Parser::ParseArg(val, args...); }) {}*/
	/*OptionDesc(const std::array<char,8> &option, const char* description, int& arg)
		: Option(option), Description(description), Handler([arg](const char* val) { Parser::ParseArg(val, arg); }) {}*/
	#define MkOpt(name, desc, args...) { name, desc, [](const OptionDesc& opt, const char* value) { opt.ParseArg(value, args); } }
	void ParseArg(const char* value, unsigned* param) const;
	void ParseArg(const char* value, unsigned* param, unsigned min, unsigned max) const;
	void ParseArg(const char* value, bool* param) const;
	void ParseArg(const char* value, double* param) const;
	void ParseArg(const char* value, const char** param) const;
	template <typename T>
	void ParseArg(const char* value, T* param, const T constvalue) const;
};

template <typename T>
void OptionDesc::ParseArg(const char* value, T* param, const T constvalue) const
{	if (*value)
		die(42, "Option %s does not expect a parameter.", Option);
	*param = constvalue;
}

#endif // PARSER_H_
