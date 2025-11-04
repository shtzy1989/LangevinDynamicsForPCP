#ifndef MD_INPUT_PARSER_H
#define MD_INPUT_PARSER_H

#include "cppheader.h"
#include "text.h"
#include "argumentParser.h"
#include "streamer.h"

namespace INPUTPARSER{
    enum class InputOptionType{
        Whatever = -2,
        CanBeLeftOut = -1,
        MustAndOnce = 0,
        Multiple = 1
    };

    enum class InputNumberRequirement{
        Whatever = -1,
        Equal = 0,
        EqualOrGreaterThan = 1
    };

    enum class InputType{
        Input = 0,
        Output = 1,
        Parameter = 2
    };

    text InputTypeToText(InputType value);

    class InputEntryValue{
    public:
        text m_ValueString;
        std::vector<text> m_ValueVector;
    };

    class InputDefinedEntry{
    public:
        // option
        text m_Option;
        text m_HelpText;
        // default value
        InputEntryValue m_DefaultValue;
        // number requriement
        InputNumberRequirement m_InputNumberRequirement;
        int m_NumberOfValues;
        // type
        InputOptionType m_InputOptionType;
        InputType m_InputType;
    public:
        InputDefinedEntry(
            text option, 
            text helpText, 
            text defaultValue, 
            InputNumberRequirement numberRequirement = InputNumberRequirement::EqualOrGreaterThan,
            int numberOfValues = 1,
            InputOptionType optionType = InputOptionType::MustAndOnce,
            InputType inputType = InputType::Parameter);
    };

    class InputEntry{
    public:
        std::vector<InputEntryValue> m_InputEntryValue;
    };

    class InputParser{
    protected:
        std::vector<InputDefinedEntry> m_DefinedEntry;
        std::vector<bool> m_EntrySpecified;
        text m_InputFileName;
        std::map<text, text> m_VariableList;
        std::map<text, text> m_DefineList;
        std::vector<InputEntry> m_InputEntry;
        std::map<text, int> m_OptionMap;
    public:
        void AddOption(text option, 
            text helpText, 
            text defaultValue, 
            InputNumberRequirement numberRequirement = InputNumberRequirement::EqualOrGreaterThan,
            int numberOfValues = 1,
            InputOptionType optionType = InputOptionType::MustAndOnce,
            InputType inputType = InputType::Parameter);
        bool Parse(int *argc, char **argv[], bool bPrint = true);
        size_t GetNumberOfOption();
        std::vector<text> GetOptionValue(text name, int valueID = 0, bool bPrint = true); // return each option value as an array
        std::vector<text> GetOptionValueListAsString(text name); // return the indexID'th value of options as an array
        std::vector<std::vector<text> > GetAllOptionValue(text name);
        std::vector<text> GetOptionValueListAsText(text name, int indexID); // return the indexID'th value of options as an array
        std::vector<double> GetOptionValueListAsDouble(text name, int indexID); // return the indexID'th value of options as an array
        std::vector<int> GetOptionValueListAsInt(text name, int indexID); // return the indexID'th value of options as an array
        std::vector<float> GetOptionValueListAsFloat(text name, int indexID); // return the indexID'th value of options as an array
        std::vector<long> GetOptionValueListAsLong(text name, int indexID); // return the indexID'th value of options as an array
        std::vector<size_t> GetOptionValueListAsSizeT(text name, int indexID); // return the indexID'th value of options as an array
        text GetOptionValueString(text name, int valueID = 0, bool bPrint = true);
        size_t GetNumberOfOptionValue(text name, bool bPrint = true);
        bool GetOptionSpecified(text name, bool bPrint = true);
        void PrintAllValues(FILE* fout);
    protected:
        bool ReadInputFile(text filename, bool bPrint, bool bOutputMode);
        bool PrintEmptyInputFile(text filename, bool bPrint);
        bool CheckDefinedEntry(bool bPrint);
        // parse INCLUDE, FOR, \, DEFINE
        bool OpenInput(text sinput, std::vector<text>* content, bool bPrint);
        std::tuple<bool, text, text> IsLineAssignment(text line);
        void GenerateMultipleLoopContent(
            std::vector<int>* startList,
            std::vector<int>* endList,
            std::vector<int>* incrementList,
            std::vector<text>* formatList,
            std::vector<text>* symbolList,
            int iloop,
            int nloop,
            text subline,
            std::vector<text>* content
            );
        bool ProcessInputLines(std::vector<text>& content, bool bPrint, bool bOutputMode);
        bool AddInputValue(int optionIndex, text value, bool bPrint);
    };
};

#endif
