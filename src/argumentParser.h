#ifndef MD_ARGUMENT_PARSER_H
#define MD_ARGUMENT_PARSER_H

#include "cppheader.h"
#include "text.h"

enum ArgumentAction{
    OnceOnly,
    Overwrite,
    Push,
};

class ArgumentEntry{
public:
    bool m_bPositional;
    int m_NumberOfValues;
    text m_OptionShort;
    text m_OptionLong;
    text m_HelpText;
    std::vector<text> m_DefaultValue;
    ArgumentAction m_ArgumentAction;
    bool m_bRequired;
    // results
    std::vector<std::vector<text> > m_Value;
    bool m_Specified;
public:
    ArgumentEntry(bool positional, int numberOfValues, text optionShort, text optionLong, text helpText, std::vector<text> defaultValue, bool bRequired, ArgumentAction argumentAction);
};

class ArgumentParser{
protected:
    std::vector<ArgumentEntry> m_DefinedEntry;
    std::map<text, size_t> m_OptionalEntryMap;
    std::map<text, size_t> m_TotalEntryMap;
    std::vector<int> m_PositionalEntryList;
    std::vector<int> m_OptionalEntryList;
public:
    void AddPositionalArgument(text optionShort, text helpText);
    void AddOptionalArgument(text optionShort, text optionLong, text helpText, int numberofValues, std::vector<text> defaultValue, bool bRequired = false, ArgumentAction ArgumentAction = ArgumentAction::OnceOnly);
    bool Parse(int *argc, char **argv[]);
    size_t GetNumberOfPositionalArgument();
    size_t GetNumberOfOptionalArgument();
    std::vector<text> GetOptionValue(text name);
    std::vector<text> GetOptionValue(text name, int valueID);
    size_t GetNumberOfOptionValue(text name);
    bool GetOptionSpecified(text name);
    void PrintAllValues();
    void PrintHelp();
private:
    bool Parse0(int *argc, char **argv[]);
};

#endif