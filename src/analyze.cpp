#include "analyze.h"

#include "parser.h"
#include "utils.h"

#include <cstdio>
#include <cstring>
#include <climits>
#include <cstdarg>
#include <cerrno>
#include <vector>
#include <algorithm>
#include <thread>

using namespace std;


/// current configuration
static Config Cfg;

/// Table of configuration parameters - MUST BE ORDERED BY NAME!
static constexpr const reference<const OptionDesc> OptMap[] =
{	MkDOp("ainc",   "incremental mode", Cfg.incremental, true)
,	MkOpt("al",     "average over multiple cycles of samples", Cfg.addloop)
,	MkSet("aref",   "append reference data file (instead of overwriting)", Cfg.refmode, "a")
,	MkOpt("bin",    "average FFT channels", Cfg.binsz)
,	MkOpt("bn",     nullptr, Cfg.N) // for compatibility
,	MkDOp("dcfg",   "dump effective configuration to file", Cfg.cfgout, "current.cfg")
,	MkDOp("diff",   "differential mode, i.e. I(t) = ch.2 - ch.1", Cfg.diffmode, true)
,	MkOpt("exec",   nullptr, Cfg.plot.shell)
,	MkOpt("famax",  "upper frequency range for LCR analysis", Cfg.famax)
,	MkOpt("famin",  "lower frequency range for LCR analysis", Cfg.famin)
,	MkOpt("fbin",   "average FFT channels with logarithmic bandwidth", Cfg.fbinsc)
,	MkSet("ff32",   "32 bit floating point format", Cfg.floatsamp, true)
,	MkOpt("fftlen", "analysis block size", Cfg.N)
,	MkSet("fi16",   "16 bit integer format", Cfg.floatsamp, false)
,	MkOpt("finc",   "linear increment for used FFT channels", Cfg.f_inc, 1., +std::numeric_limits<double>::infinity())
,	MkOpt("flog",   "logarithmic increment for used FFT channels", Cfg.f_log, 1., +std::numeric_limits<double>::infinity())
,	MkOpt("fmax",   "upper frequency range for analysis", Cfg.fmax, 1E-16, 1E16)
,	MkOpt("fmin",   "lower frequency range for analysis", Cfg.fmin, 1E-16, 1E16)
,	MkOpt("fsamp",  "sampling frequency, 48k by default", Cfg.srate)
,	MkOpt("gain",   "output gain in dB FSR", Cfg.outgain)
,	MkDOp("gg",     "generate gain calibration file", Cfg.gainoutfile, "gain.dat")
,	MkDOp("gr",     "use gain calibration file", Cfg.gaininfile, "gain.dat")
,	MkSet("h/f",    "use 1/f weight", Cfg.weightfn, &Config::Get1_fWeight)
,	MkOpt("harm",   "take harmonics into account", Cfg.harmonic, 1U, HA_MAX)
,	MkSet("hd",     "use weight function for differential input mode", Cfg.weightfn, &Config::GetWeightD)
,	MkSet("he",     "disable weight function", Cfg.weightfn, &Config::GetConstWeight)
,	MkSet("help",   "show this help", Cfg.help, true)
,	MkDOp("in",     "name of input file to analyze, stdin by default", Cfg.infile, "-")
,	MkOpt("initcmd","execute shell command before any data processing", Cfg.init.shell)
,	MkOpt("initout","print string to stdout before any data processing", Cfg.init.out)
,	MkOpt("ln",     nullptr, Cfg.loops)
,	MkSet("loop",   "infinite number of loops", Cfg.loops, 0U)
,	MkOpt("loops",  "number of loops", Cfg.loops)
,	MkOpt("lp",     "pause at matrix calibration", Cfg.lpause, 0., 100.)
,	MkSet("mfft",   "enable operation mode FFT", Cfg.mfft, true)
,	MkDOp("minmax", "Show minimum and maxoimum values of input data", Cfg.minmax, true)
,	MkSet("mpca",   "enable operation mode PCA", Cfg.mpca, true)
,	MkDOp("mst",    nullptr, Cfg.stereo, true)
,	MkSet("msweep", "use sweep mode", Cfg.sweep, true)
,	MkSet("mxy",    "enable operation mode XY (preliminary)", Cfg.mxy, true)
,	MkDOp("nohdr",  "do not generate WAV header", Cfg.nohdr, true)
,	MkOpt("olc",    "column to overwrite numerator", Cfg.overwrt[0].column)
,	MkOpt("olf",    "file name to overwrite numerator", Cfg.overwrt[0].file)
,	MkOpt("orc",    "column to overwrite denominator (reference)", Cfg.overwrt[1].column)
,	MkOpt("orf",    "file name to overwrite denominator (reference)", Cfg.overwrt[1].file)
,	MkDOp("out",    "name of reference output file, stdout by default", Cfg.outfile, "-")
,	MkOpt("pch",    "purge first frequency channels", Cfg.purgech)
,	MkDOp("phcc",   "fit group delay", Cfg.crosscorr, true)
,	MkOpt("phl",    "subtract constant group delay", Cfg.linphase)
,	MkOpt("plot",   nullptr, Cfg.plot.out)
,	MkOpt("plotcmd","execute shell command after data available", Cfg.plot.shell)
,	MkOpt("plotout","print string to stdout after data available", Cfg.plot.out)
,	MkOpt("postcmd","execute shell command after completion", Cfg.post.shell)
,	MkOpt("postout","print string to stdout after completion", Cfg.post.out)
,	MkOpt("predelay","setup delay in FFT cycles", Cfg.predelay, 0., 100.)
,	MkOpt("psa",    "discard first samples", Cfg.discsamp)
,	MkDOp("pte",    "read input data till the end", Cfg.disctrail, true)
,	MkOpt("rref",   "reference resistor", Cfg.rref)
,	MkDOp("rspec",  "read frequency domain reference data from file", Cfg.rspecfile, "spectrum.dat")
,	MkOpt("scale",  "noise type", Cfg.scalepow)
,	MkOpt("setupcmd","execute shell command at program start", Cfg.setup.shell)
,	MkOpt("setupout","print string to stdout at program start", Cfg.setup.out)
,	MkDOp("stereo", "two channel mode", Cfg.stereo, true)
,	MkDOp("symmout","symmetric output", Cfg.symmout, true)
,	MkDOp("sync",   "synchronize cycles before start", Cfg.sync, 2U)
,	MkOpt("synclvl","minimum correlation ratio of successful synchronization", Cfg.synclevel, 0., 1.)
,	MkDOp("wd",     "(over)write FFT data file on the fly", Cfg.datafile, "data.dat")
,	MkOpt("win",    "select window function [0..5]", Cfg.winfn, 0U, 5U)
,	MkDOp("wraw",   "write raw input data to file", Cfg.rawfile, "raw.dat")
,	MkDOp("wref",   "write time domain reference data to file", Cfg.reffile, "ref.dat")
,	MkDOp("wspec",  "write frequency domain reference data to file", Cfg.specfile, "spectrum.dat")
,	MkDOp("wsrc",   "write source data file", Cfg.srcfile, "source.dat")
,	MkDOp("wwin",   "write window function to file", Cfg.windowfile, "window.dat")
,	MkDOp("xch",    "swap input channels L<->R", Cfg.swapch, true)
,	MkDOp("xb" ,    "swap bytes", Cfg.swapbytes, true)
,	MkDOp("zg",     "generate matrix calibration file", Cfg.zerooutfile, "zero.dat")
,	MkDOp("zn",     "normalize amplitudes", Cfg.normalize, true)
,	MkDOp("zr",     "use matrix calibration file", Cfg.zeroinfile, "zero.dat")
};


void action::execute() const
{	::execute(shell);
	if (out)
	{	puts(out);
		fflush(stdout);
	}
}


static void do_parallel(const vector<ITask*>& tasks)
{	switch (tasks.size())
	{case 1:
		tasks[0]->Run();
	 case 0:
		break;
	 default:
		vector<thread> threads;
		threads.reserve(tasks.size() - 1);
		for (auto task : tasks)
			if (threads.size() < tasks.size() -1)
				threads.emplace_back(&ITask::Run, task);
			else
				task->Run();
		for (auto& thread : threads)
			thread.join();
	}
}


int main(int argc, char* argv[])
{
	// parse cmdl
	{	Parser parser(OptMap);
		while (--argc)
			parser.HandleArg(*++argv);

		if (Cfg.help || (!Cfg.infile &&!Cfg.outfile && !Cfg.windowfile && !Cfg.reffile && !Cfg.specfile))
		{	fputs("usage: analyze <options>  -  see documentation for more details.\n", stderr);
			parser.PrintHelp();
			return 48;
		}

		// Verify config
		if (Cfg.N > N_MAX)
			die(32, "FFT Length too large.");
		if (Cfg.N & (Cfg.N - 1))
			fputs("Warning: FFT size should be a power of to to achieve optimum performance.\n", stderr);
		if (Cfg.overwrt[0].file && Cfg.overwrt[1].file)
			Cfg.infile = NULL;
		if (Cfg.mxy & (Cfg.mpca|Cfg.mfft))
			die(34, "Invalid combination of measurement modes, e.g. FFT and XY.");
		if (Cfg.outfile && strcmp(Cfg.outfile,"-") == 0 && (Cfg.init.out || Cfg.plot.out || Cfg.post.out))
			die(34, "Options initout/plotout/postout cannot be used when writing PCM data to stdout.");
		if (!Cfg.loops && (Cfg.sweep || Cfg.zerooutfile))
			die(34, "Cannot use infinite loop in sweep or zero calibration mode.");

		if (Cfg.fmin > Cfg.fmax)
			die(34, "Frequency range [%g,%g) empty.", Cfg.fmin, Cfg.fmax);
		if (Cfg.famin > Cfg.famax)
			die(34, "Analysis frequency range [%g,%g) empty.", Cfg.famin, Cfg.famax);
		if (Cfg.fmin < (double)Cfg.srate/Cfg.N)
			Cfg.fmin = (double)Cfg.srate/Cfg.N;
		if (Cfg.fmax > Cfg.srate/2)
			Cfg.fmax = Cfg.srate/2;
		if (Cfg.famin < Cfg.fmin)
			Cfg.famin = Cfg.fmin;
		if (Cfg.famax > Cfg.fmax)
			Cfg.famax = Cfg.fmax;

		if (Cfg.cfgout)
			parser.WriteConfig(FILEguard(Cfg.cfgout, "wt"));
	}

	Cfg.setup.execute();

	srand(clock());

	vector<ITask*> tasks;
	tasks.reserve(2);
	unique_ptr<AnalyzeOut> source(new AnalyzeOut(Cfg));
	if (source->Setup())
		tasks.emplace_back(source.get());
	unique_ptr<AnalyzeIn> drain(new AnalyzeIn(Cfg, *source));
	if (drain->Setup())
		tasks.emplace_back(drain.get());

	Cfg.init.execute();

	do_parallel(tasks);

	return 0;
}
