#ifndef PARSER_H_
#define PARSER_H_

#include "utils.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

// get rid of many false positives with non-virtual dtor warning.
//#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"

/// Abstract base class of command line argument descriptor.
struct OptionDesc
{	const cstring<8> Option;
	const char* const Description;
	constexpr    OptionDesc(cstring<8> opt, const char* desc) : Option(opt), Description(desc) {}
	virtual void Parse(const char* val) const = 0;
	virtual void Print() const = 0;
	virtual void Write(FILE* out) const = 0;

	// implementing classes
	struct SetBase;
	template <typename T>	struct Set;
	template <typename T>	struct Opt{}; // invalid type, see specializations below
	template <typename T>	struct DefOpt;
	template <typename T>	struct LimitOpt{}; // invalid type, see specializations below

	struct Comparer
	{	bool operator()(const OptionDesc& opt, const char* val) { return std::strncmp(opt.Option.data(), val, 8) < 0; }
		bool operator()(const char* val, const OptionDesc& opt) { return std::strncmp(val, opt.Option.data(), 8) < 0; }
	};
 protected:
	void PrintPreamble() const;
	virtual size_t PrintParam() const = 0;
};

struct OptionDesc::SetBase : OptionDesc
{protected:
	constexpr SetBase(cstring<8> opt, const char* desc) : OptionDesc(opt, desc) {}
	virtual size_t PrintParam() const;
	void CheckNoVal(const char* val) const;
	void DoPrint(bool active) const;
	void DoWrite(FILE* out) const;
};
template <typename T>
struct OptionDesc::Set final : OptionDesc::SetBase
{	T& Param;
	const T Value;
	constexpr Set(cstring<8> opt, const char* desc, T& par, const T& val) : SetBase(opt, desc), Param(par), Value(val) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
	virtual void Write(FILE* out) const;
};
template <typename T>
void OptionDesc::Set<T>::Parse(const char* val) const
{	CheckNoVal(val);
	Param = Value;
}
template <typename T>
void OptionDesc::Set<T>::Print() const
{	DoPrint(Param == Value);
}
template <typename T>
void OptionDesc::Set<T>::Write(FILE* out) const
{	if (Param == Value)
		DoWrite(out);
}

template <>
struct OptionDesc::Opt<bool> : OptionDesc
{	bool& Param;
	constexpr Opt(cstring<8> opt, const char* desc, bool& par) : OptionDesc(opt, desc), Param(par) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
	virtual void Write(FILE* out) const;
 protected:
	virtual size_t PrintParam() const;
	static void PrintValue(bool value);
};
template <>
struct OptionDesc::Opt<unsigned> : OptionDesc
{	unsigned& Param;
	constexpr Opt(cstring<8> opt, const char* desc, unsigned& par) : OptionDesc(opt, desc), Param(par) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
	virtual void Write(FILE* out) const;
 protected:
	virtual size_t PrintParam() const;
	static void PrintValue(unsigned value);
};
template <>
struct OptionDesc::Opt<double> : OptionDesc
{	double& Param;
	constexpr Opt(cstring<8> opt, const char* desc, double& par) : OptionDesc(opt, desc), Param(par) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
	virtual void Write(FILE* out) const;
 protected:
	virtual size_t PrintParam() const;
	static void PrintValue(double value);
};
template <>
struct OptionDesc::Opt<const char*> : OptionDesc
{	const char*& Param;
	constexpr Opt(cstring<8> opt, const char* desc, const char*& par) : OptionDesc(opt, desc), Param(par) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
	virtual void Write(FILE* out) const;
 protected:
	virtual size_t PrintParam() const;
	static void PrintValue(const char* value);
};

template <typename T>
struct OptionDesc::DefOpt final : OptionDesc::Opt<T>
{	const T Default;
	constexpr DefOpt(cstring<8> opt, const char* desc, T& par, T def) : OptionDesc::Opt<T>(opt, desc, par), Default(def) {}
	virtual void Parse(const char* val) const;
	virtual void Print() const;
};

template <>
struct OptionDesc::LimitOpt<unsigned> final : OptionDesc::Opt<unsigned>
{	const unsigned Min;
	const unsigned Max;
	constexpr LimitOpt(cstring<8> opt, const char* desc, unsigned& par, unsigned min, unsigned max) : OptionDesc::Opt<unsigned>(opt, desc, par), Min(min), Max(max) {}
	virtual void Parse(const char* val) const;
};
template <>
struct OptionDesc::LimitOpt<double> final : OptionDesc::Opt<double>
{	const double Min;
	const double Max;
	constexpr LimitOpt(cstring<8> opt, const char* desc, double& par, double min, double max) : OptionDesc::Opt<double>(opt, desc, par), Min(min), Max(max) {}
	virtual void Parse(const char* val) const;
};

template <typename T>
constexpr OptionDesc::Set<T> MkSet(cstring<8> opt, const char* desc, T& par, T set)
{	return OptionDesc::Set<T>(opt, desc, par, set);
}
template <typename T>
constexpr OptionDesc::Opt<T> MkOpt(cstring<8> opt, const char* desc, T& par)
{	return OptionDesc::Opt<T>(opt, desc, par);
}
template <typename T>
constexpr OptionDesc::LimitOpt<T> MkOpt(cstring<8> opt, const char* desc, T& par, T min, T max)
{	return OptionDesc::LimitOpt<T>(opt, desc, par, min, max);
}
template <typename T>
constexpr OptionDesc::DefOpt<T> MkDOp(cstring<8> opt, const char* desc, T& par, T def)
{	return OptionDesc::DefOpt<T>(opt, desc, par, def);
}

template <typename T>
void OptionDesc::DefOpt<T>::Parse(const char* val) const
{	if (!val)
		Opt<T>::Param = Default;
	else
		Opt<T>::Parse(val);
}
template <typename T>
void OptionDesc::DefOpt<T>::Print() const
{	Opt<T>::PrintPreamble();
	fputs(" (current: ", stderr);
	Opt<T>::PrintValue(Opt<T>::Param);
	fputs(", default: ", stderr);
	Opt<T>::PrintValue(Default);
	fputc(')', stderr);
}

class Parser
{	const reference<const OptionDesc>* const ArgMap;
	const reference<const OptionDesc>* const ArgMapEnd;
 public:
	template <size_t N>
	/// Create option parser for list of options
	/// @param argmap List of \see OptionDesc \b ordered by name.
	Parser(const reference<const OptionDesc> (&argmap)[N]) : ArgMap(argmap), ArgMapEnd(argmap + N) {}
	/// Apply option argument to the configuration.
	void HandleArg(const char* arg);
	/// Print help about all (documented) options to \c stderr.
	void PrintHelp() const;
	///Write all (documented) configuration options to a file.
	void WriteConfig(FILE* out) const;
};

#endif // PARSER_H_
