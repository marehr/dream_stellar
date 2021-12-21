#include "stellar.arg_parser.hpp"

#include <seqan/seq_io.h>

namespace stellar
{

namespace app
{

///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
ArgumentParser::ParseResult
_parseOptions(ArgumentParser const & parser, StellarOptions & options)
{
    getArgumentValue(options.databaseFile, parser, 0);
    getArgumentValue(options.queryFile, parser, 1);

    // output options
    getOptionValue(options.outputFile, parser, "out");
    getOptionValue(options.disabledQueriesFile, parser, "outDisabled");
    getOptionValue(options.noRT, parser, "no-rt");

    CharString tmp = options.outputFile;
    toLower(tmp);

    if (endsWith(tmp, ".gff"))
        options.outputFormat = "gff";
    else if (endsWith(tmp, ".txt"))
        options.outputFormat = "txt";

    // main options
    getOptionValue(options.qGram, parser, "kmer");
    getOptionValue(options.minLength, parser, "minLength");
    getOptionValue(options.epsilon, parser, "epsilon");
    getOptionValue(options.xDrop, parser, "xDrop");
    getOptionValue(options.alphabet, parser, "alphabet");
    getOptionValue(options.threadCount, parser, "threads");

    if (isSet(parser, "forward") && !isSet(parser, "reverse"))
        options.reverse = false;
    if (!isSet(parser, "forward") && isSet(parser, "reverse"))
        options.forward = false;

    getOptionValue(options.fastOption, parser, "verification");
    getOptionValue(options.disableThresh, parser, "disableThresh");
    getOptionValue(options.numMatches, parser, "numMatches");
    getOptionValue(options.compactThresh, parser, "sortThresh");
    getOptionValue(options.maxRepeatPeriod, parser, "repeatPeriod");
    getOptionValue(options.minRepeatLength, parser, "repeatLength");
    getOptionValue(options.qgramAbundanceCut, parser, "abundanceCut");

    getOptionValue(options.verbose, parser, "verbose");

    if (isSet(parser, "kmer") && options.qGram >= 1 / options.epsilon)
    {
        std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    if (options.numMatches > options.compactThresh)
    {
        std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

///////////////////////////////////////////////////////////////////////////////
// Set-Up of Argument Parser
void _setParser(ArgumentParser & parser)
{
    setShortDescription(parser, "the SwifT Exact LocaL AligneR");
    setDate(parser, SEQAN_DATE);
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setCategory(parser, "Local Alignment");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIFASTA FILE 1\\fP> <\\fIFASTA FILE 2\\fP>");

    addDescription(parser,
                   "STELLAR implements the SWIFT filter algorithm (Rasmussen et al., 2006) "
                   "and a verification step for the SWIFT hits that applies local alignment, "
                   "gapped X-drop extension, and extraction of the longest epsilon-match.");
    addDescription(parser,
                   "Input to STELLAR are two files, each containing one or more sequences "
                   "in FASTA format. Each sequence from file 1 will be compared to each "
                   "sequence in file 2. The sequences from file 1 are used as database, the "
                   "sequences from file 2 as queries.");
    addDescription(parser, "(c) 2010-2012 by Birte Kehr");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 1"));
    setValidValues(parser, 0, "fa fasta");  // allow only fasta files as input
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 2"));
    setValidValues(parser, 1, "fa fasta");  // allow only fasta files as input

    // Add threads option.
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setDefaultValue(parser, "threads", "1");

    addSection(parser, "Main Options");

    addOption(parser, ArgParseOption("e", "epsilon", "Maximal error rate (max 0.25).", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e", "0.05");
    setMinValue(parser, "e", "0.0000001");
    setMaxValue(parser, "e", "0.25");
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length of epsilon-matches.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "l", "100");
    setMinValue(parser, "l", "0");
    addOption(parser, ArgParseOption("f", "forward", "Search only in forward strand of database."));
    addOption(parser, ArgParseOption("r", "reverse", "Search only in reverse complement of database."));
    addOption(parser, ArgParseOption("a", "alphabet",
                                     "Alphabet type of input sequences (dna, rna, dna5, rna5, protein, char).",
                                     ArgParseArgument::STRING));
    setValidValues(parser, "a", "dna dna5 rna rna5 protein char");
    addOption(parser, ArgParseOption("v", "verbose", "Set verbosity mode."));

    addSection(parser, "Filtering Options");

    addOption(parser, ArgParseOption("k", "kmer", "Length of the q-grams (max 32).", ArgParseArgument::INTEGER));
    setMinValue(parser, "k", "1");
    setMaxValue(parser, "k", "32");
    addOption(parser, ArgParseOption("rp", "repeatPeriod",
                                     "Maximal period of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rp", "1");
    addOption(parser, ArgParseOption("rl", "repeatLength",
                                     "Minimal length of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rl", "1000");
    addOption(parser, ArgParseOption("c", "abundanceCut", "k-mer overabundance cut ratio.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "c", "1");
    setMinValue(parser, "c", "0");
    setMaxValue(parser, "c", "1");

    addSection(parser, "Verification Options");

    addOption(parser, ArgParseOption("x", "xDrop", "Maximal x-drop for extension.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "x", "5");
    addOption(parser, ArgParseOption("vs", "verification", "Verification strategy: exact or bestLocal or bandedGlobal",
                                     ArgParseArgument::STRING));
    //addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
    //addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
    //addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
    setDefaultValue(parser, "vs", "exact");
    setValidValues(parser, "vs", "exact bestLocal bandedGlobal");
    addOption(parser, ArgParseOption("dt", "disableThresh",
                                     "Maximal number of verified matches before disabling verification for one query "
                                     "sequence (default infinity).", ArgParseArgument::INTEGER));
    setMinValue(parser, "dt", "0");
    addOption(parser, ArgParseOption("n", "numMatches",
                                     "Maximal number of kept matches per query and database. If STELLAR finds more matches, "
                                     "only the longest ones are kept.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "n", "50");
    addOption(parser, ArgParseOption("s", "sortThresh",
                                     "Number of matches triggering removal of duplicates. Choose a smaller value for saving "
                                     "space.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "s", "500");

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "out", "Name of output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "o", "gff txt");
    setDefaultValue(parser, "o", "stellar.gff");
    addOption(parser, ArgParseOption("od", "outDisabled",
                                     "Name of output file for disabled query sequences.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "outDisabled", seqan::SeqFileOut::getFileExtensions());
    setDefaultValue(parser, "od", "stellar.disabled.fasta");
    addOption(parser, ArgParseOption("no-rt", "suppress-runtime-printing", "Suppress printing running time."));
    hideOption(parser, "no-rt");

    addTextSection(parser, "References");
    addText(parser, "Kehr, B., Weese, D., Reinert, K.: STELLAR: fast and exact local alignments. BMC Bioinformatics, "
                    "12(Suppl 9):S15, 2011.");
}

} // namespace stellar::app

} // namespace stellar
