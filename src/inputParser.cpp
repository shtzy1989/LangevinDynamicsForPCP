#include "inputParser.h"
#include "MOLUTILITY.h"

namespace INPUTPARSER{
    text InputTypeToText(InputType value){
        switch(value){
            case InputType::Input:     return "Input";
            case InputType::Output:    return "Output";
            case InputType::Parameter: return "Parameter";
            default: return "";
        }
    };

    InputDefinedEntry::InputDefinedEntry(
        text option, 
        text helpText, 
        text defaultValue, 
        InputNumberRequirement numberRequirement,
        int numberOfValues,
        InputOptionType optionType,
        InputType inputType){
        // code
        m_Option = option;
        m_HelpText = helpText;
        m_DefaultValue.m_ValueString = defaultValue;
        m_DefaultValue.m_ValueVector = m_DefaultValue.m_ValueString.split();
        m_InputNumberRequirement = numberRequirement;
        m_NumberOfValues = numberOfValues;
        m_InputOptionType = optionType;
        m_InputType = inputType;
    }

    void InputParser::AddOption(text option, 
        text helpText, 
        text defaultValue, 
        InputNumberRequirement numberRequirement,
        int numberOfValues,
        InputOptionType optionType,
        InputType inputType){
        // code
        m_DefinedEntry.emplace_back(
            option, helpText, defaultValue, numberRequirement, 
            numberOfValues, optionType, inputType);
    }
    bool InputParser::Parse(int *argc, char **argv[], bool bPrint){
        ArgumentParser argumentParser;
        argumentParser.AddOptionalArgument("-i", "--input", "Input file name", 1, {"none"}, false);
        argumentParser.AddOptionalArgument("-o", "--output", "Input file name to be generated", 1, {"none"}, false);
        argumentParser.AddOptionalArgument("-p", "--print", "Print input file content", 0, {}, false);
        argumentParser.Parse(argc, argv);

        if( !CheckDefinedEntry(bPrint) ){
            // exit(0);
            return false;
        }

        m_InputFileName = argumentParser.GetOptionValue("-i")[0];
        text inputFileNameToBeGenerated = argumentParser.GetOptionValue("-o")[0];
        if( argumentParser.GetOptionSpecified("-i") && 
            argumentParser.GetOptionSpecified("-o") ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): only one of -i and -o can be specified\n");
            // exit(0);
            return false;
        }else if( argumentParser.GetOptionSpecified("-i") ){
            bool bOK = ReadInputFile(m_InputFileName, bPrint, false);
            if( bOK ){
                if( argumentParser.GetOptionSpecified("-p") ){
                    if( bPrint ) PrintAllValues(stderr);
                    // exit(0);
                    return false;
                }else{
                    return true;
                }
            }else{
                return false;
            }
        }else if( argumentParser.GetOptionSpecified("-o") ){
            bool bOK = PrintEmptyInputFile(inputFileNameToBeGenerated, bPrint);
            // exit(0);
            return false;
        }else{
            if( bPrint ){
                argumentParser.PrintHelp();
            }
            // exit(0);
            return false;
        }
    };
    size_t InputParser::GetNumberOfOption(){
        return m_DefinedEntry.size();
    }
    std::vector<text> InputParser::GetOptionValue(text name, int valueID, bool bPrint){
        int optionIndex = m_OptionMap[name] - 1;
        if( optionIndex == -1 ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): undefined option (%s)\n", name.c_str());
            // exit(0);
        }else{
            return m_InputEntry[optionIndex].m_InputEntryValue[valueID].m_ValueVector;
        }
        return {};
    }
    std::vector<text> InputParser::GetOptionValueListAsString(text name){
        std::vector<text> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValueString(name, i));
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<std::vector<text> > InputParser::GetAllOptionValue(text name){
        int numberOfOptionValues = GetNumberOfOptionValue(name);
        std::vector<std::vector<text> > result(numberOfOptionValues);
        for(int i=0;i<numberOfOptionValues;i++){
            result[i] = GetOptionValue(name, i);
        }
        return result;
    }
    std::vector<text> InputParser::GetOptionValueListAsText(text name, int indexID){
        std::vector<text> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID]);
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<double> InputParser::GetOptionValueListAsDouble(text name, int indexID){
        std::vector<double> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID].to_double());
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<int> InputParser::GetOptionValueListAsInt(text name, int indexID){
        std::vector<int> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID].to_int());
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<float> InputParser::GetOptionValueListAsFloat(text name, int indexID){
        std::vector<float> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID].to_float());
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<long> InputParser::GetOptionValueListAsLong(text name, int indexID){
        std::vector<long> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID].to_long());
        }
        return result;
    }; // return the indexID'th value of options as an array
    std::vector<size_t> InputParser::GetOptionValueListAsSizeT(text name, int indexID){
        std::vector<size_t> result;
        for(int i=0;i<GetNumberOfOptionValue(name);i++){
            result.push_back(GetOptionValue(name, i)[indexID].to_unsigned_long());
        }
        return result;
    }; // return the indexID'th value of options as an array
    text InputParser::GetOptionValueString(text name, int valueID, bool bPrint){
        int optionIndex = m_OptionMap[name] - 1;
        if( optionIndex == -1 ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): undefined option (%s)\n", name.c_str());
            // exit(0);
        }else{
            return m_InputEntry[optionIndex].m_InputEntryValue[valueID].m_ValueString;
        }
        return "";
    }
    size_t InputParser::GetNumberOfOptionValue(text name, bool bPrint){
        int optionIndex = m_OptionMap[name] - 1;
        if( optionIndex == -1 ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): undefined option (%s)\n", name.c_str());
            // exit(0);
        }else{
            return m_InputEntry[optionIndex].m_InputEntryValue.size();
        }
        return 0;
    }
    bool InputParser::GetOptionSpecified(text name, bool bPrint){
        int optionIndex = m_OptionMap[name] - 1;
        if( optionIndex == -1 ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): undefined option (%s)\n", name.c_str());
            // exit(0);
        }else{
            return m_EntrySpecified[optionIndex];
        }
        return false;
    }
    void InputParser::PrintAllValues(FILE* fout){
        for(int i=0;i<m_InputEntry.size();i++){
            for(int j=0;j<m_InputEntry[i].m_InputEntryValue.size();j++){
                fprintf(stderr, "%-20s %3d %s\n", 
                    m_DefinedEntry[i].m_Option.c_str(), j, 
                    m_InputEntry[i].m_InputEntryValue[j].m_ValueString.c_str());
            }
        }
    }
    bool InputParser::ReadInputFile(text filename, bool bPrint, bool bOutputMode){
        std::vector<text> content;
        if( !OpenInput(filename, &content, bPrint) ){
            // if( bExitOnFail ) exit(0);
            return false;
        }else{
            // process
            bool bOK = ProcessInputLines(content, bPrint, bOutputMode);
            if( !bOK ){
                // if( bExitOnFail ) exit(0);
            }
            return bOK;
        }
    }
    bool InputParser::PrintEmptyInputFile(text filename, bool bPrint){
        bool bUseExistInputValue = false;
        if( MOLUTILITY::file_exist(filename.string()) ){
            // read existing input file
            bUseExistInputValue = ReadInputFile(filename, false, true);
        }
        TextFileWriter writer;
        if( !writer.Open(filename, true) ){
            fprintf(stderr, "Error (InputParser): Can't open empty input file (%s)\n", filename.c_str());
            return false;
        }else{
            for(int i=0;i<m_DefinedEntry.size();i++){
                writer.Write("# [%s] %s\n", 
                    InputTypeToText(m_DefinedEntry[i].m_InputType).c_str(),
                    m_DefinedEntry[i].m_HelpText.c_str());
                if( bUseExistInputValue ){
                    if( m_InputEntry[i].m_InputEntryValue.size() != 0 ){
                        for(int j=0;j<m_InputEntry[i].m_InputEntryValue.size();j++){
                            writer.Write("%-20s %s\n", 
                                m_DefinedEntry[i].m_Option.c_str(), 
                                m_InputEntry[i].m_InputEntryValue[j].m_ValueString.c_str());
                        }
                    }else{
                        writer.Write("%-20s %s\n", 
                            m_DefinedEntry[i].m_Option.c_str(), 
                            m_DefinedEntry[i].m_DefaultValue.m_ValueString.c_str());
                    }
                }else{
                    writer.Write("%-20s %s\n", 
                        m_DefinedEntry[i].m_Option.c_str(), 
                        m_DefinedEntry[i].m_DefaultValue.m_ValueString.c_str());
                }
            }
            writer.Close();
            return true;
        }
    }
    bool InputParser::CheckDefinedEntry(bool bPrint){
        bool bOK = true;
        for(int i=0;i<m_DefinedEntry.size();i++){
            if( m_DefinedEntry[i].m_InputNumberRequirement == InputNumberRequirement::Equal ){
                if( m_DefinedEntry[i].m_DefaultValue.m_ValueVector.size() != m_DefinedEntry[i].m_NumberOfValues ){
                    if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) requires %d values but got %zd default values\n",
                        m_DefinedEntry[i].m_Option.c_str(), m_DefinedEntry[i].m_NumberOfValues, m_DefinedEntry[i].m_DefaultValue.m_ValueVector.size());
                    bOK = false;
                }
            }else if( m_DefinedEntry[i].m_InputNumberRequirement == InputNumberRequirement::EqualOrGreaterThan ){
                if( m_DefinedEntry[i].m_DefaultValue.m_ValueVector.size() < m_DefinedEntry[i].m_NumberOfValues ){
                    if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) requires %d values but got %zd default values\n",
                        m_DefinedEntry[i].m_Option.c_str(), m_DefinedEntry[i].m_NumberOfValues, m_DefinedEntry[i].m_DefaultValue.m_ValueVector.size());
                    bOK = false;
                }
            }
        }
        return bOK;
    }
    // parse INCLUDE, FOR, \, DEFINE
    bool InputParser::OpenInput(text sinput, std::vector<text>* content, bool bPrint){
        bool bResult = true;
        StreamerReader reader;
        if( !reader.Open(sinput) ){
            if( bPrint ) fprintf(stderr, "Error (InputParser): can't open input file (%s)\n", 
                m_InputFileName.c_str());
            bResult = false;
        }else{
            // read lines and parse \,
            std::vector<text> lines;
            bool bAppendLine = false;
            while( reader.NextLineAvailable() ){
                text line = reader.ReadLine();
                line = line.parseRemoveComment("#");
                if( bAppendLine ){
                    lines.back() += line;
                }else{
                    lines.push_back(line);
                }
                if( line.rightEqual("\\") ){
                    lines.back().resize(lines.back().size() - 1);
                    bAppendLine = true;
                }else{
                    bAppendLine = false;
                }
            }
            for(int i=0;i<lines.size();i++){
                if( lines[i].isEmpty() ) continue;
                // apply define and variables
                for(auto j=m_DefineList.begin();j!=m_DefineList.end();j++){
                    lines[i] = lines[i].replaceAll(j->first, j->second);
                }
                for(auto j=m_VariableList.begin();j!=m_VariableList.end();j++){
                    lines[i] = lines[i].replaceAll(j->first, j->second);
                }
                auto word = lines[i].split();
                auto [bAssignment, lhs, rhs] = IsLineAssignment(lines[i]);
                if( bAssignment ){
                    m_VariableList[lhs] = rhs;
                }else if( word[0].ContentEqualCaseInsensitive("include") ){
                    if( word.size() >= 2 ){
                        for(int w=1;w<word.size();w++){
                            bResult = bResult && OpenInput(word[w], content, bPrint);
                            if( !bResult ){
                                break;
                            }
                        }
                    }else{
                        if( bPrint ) fprintf(stderr, "Error (InputParser): %s requires at least 1 argument, but now %d\n", 
                            word[0].c_str(), (int)word.size() - 1);
                        bResult = false;
                    }
                }else if( word[0].ContentEqualCaseInsensitive("define") ){
                    if( word.size() == 3 ){
                        m_DefineList[word[1]] = word[2];
                    }else{
                        if( bPrint ) fprintf(stderr, "Error (InputParser): %s requires at least 1 argument, but now %d\n", 
                            word[0].c_str(), (int)word.size() - 1);
                        bResult = false;
                    }
                }else if( word[0].ContentEqualCaseInsensitive("for") ){
                    if( word.size() >= 7 ){
                        int iStart = word[1].to_int();
                        int iEnd = word[2].to_int();
                        int iFreq = word[3].to_int();
                        text format = word[4];
                        text symbol = word[5];
                        text command = text::FromTextArray(&word, 6);
                        for(int j=iStart;j<iEnd;j+=iFreq){
                            text newLine = command.replaceAll(symbol, text::formatC(format.c_str(), j));
                            content->push_back(newLine);
                        }
                    }
                }else if( word[0].left(3).ContentEqualCaseInsensitive("for")){
                    if( word[0].substr(3, word[0].size() - 3).isInt() ){
                        int nloop = word[0].substr(3, word[0].size() - 3).to_int();
                        std::vector<int> startList(nloop);
                        std::vector<int> endList(nloop);
                        std::vector<int> incrementList(nloop);
                        std::vector<text> formatList(nloop);
                        std::vector<text> symbolList(nloop);
                        for(int i=0;i<nloop;i++){
                            startList[i] = word[i*5+1].to_int();
                            endList[i] = word[i*5+2].to_int();
                            incrementList[i] = word[i*5+3].to_int();
                            formatList[i] = word[i*5+4];
                            symbolList[i] = word[i*5+5];
                        }
                        text subline = text::FromTextArray(&word, nloop * 5 + 1);
                        GenerateMultipleLoopContent(
                            &startList,
                            &endList,
                            &incrementList,
                            &formatList,
                            &symbolList,
                            0,
                            nloop,
                            subline,
                            content);
                    }else{
                        if( bPrint ) fprintf(stderr, "Error (InputParser): for is not followed by an integer (%s)\n", 
                            word[0].substr(3, word[0].size() - 3).c_str());
                    }
                }else{
                    content->push_back(lines[i]);
                }
                if( !bResult ){
                    break;
                }
            }
            reader.Close();
        }
        return bResult;
    }
    std::tuple<bool, text, text> InputParser::IsLineAssignment(text line){
        if( line.Contain("=") ){
            auto word = line.split();
            auto wordForEqual = line.split(" \t=");
            if( wordForEqual.size() == 2 ){
                auto word0 = wordForEqual[0].split();
                auto word1 = wordForEqual[1].split();
                if( word0.size() == 1 && word1.size() == 1 ){
                    return {true, word0[0], word1[0]};
                }else{
                    return {false, "", ""};
                }
            }else{
                return {false, "", ""};
            }
        }else{
            return {false, "", ""};
        }
    }
    void InputParser::GenerateMultipleLoopContent(
        std::vector<int>* startList,
        std::vector<int>* endList,
        std::vector<int>* incrementList,
        std::vector<text>* formatList,
        std::vector<text>* symbolList,
        int iloop,
        int nloop,
        text subline,
        std::vector<text>* content
        ){
        for(int i=(*startList)[iloop];i<(*endList)[iloop];i+=(*incrementList)[iloop]){
            char buffer [256];
            text newline = subline;
            sprintf(buffer, (*formatList)[iloop].c_str(), i);
            newline = newline.replaceAll((*symbolList)[iloop], buffer);
            if( iloop == nloop - 1 ){
                content->push_back(newline);
            }else{
                GenerateMultipleLoopContent(
                    startList,
                    endList,
                    incrementList,
                    formatList,
                    symbolList,
                    iloop + 1,
                    nloop,
                    newline,
                    content);
            }
        }
    };
    bool InputParser::ProcessInputLines(std::vector<text>& content, bool bPrint, bool bOutputMode){
// for(int i=0;i<content.size();i++){
//     fprintf(stderr, "%3d %s\n", i, content[i].c_str());
// }
        for(int i=0;i<m_DefinedEntry.size();i++){
            m_OptionMap[m_DefinedEntry[i].m_Option] = i + 1;
        }

        m_InputEntry.resize(m_DefinedEntry.size());
        m_EntrySpecified.resize(m_DefinedEntry.size(), false);
        for(int i=0;i<content.size();i++){
            if( !content[i].isEmpty() ){
                auto word = content[i].split();
                int optionIndex = m_OptionMap[word[0]] - 1;
                if( optionIndex == -1 ){
                    if( bPrint ) fprintf(stderr, "Error (InputParser): undefined option (%s)\n", 
                        word[0].c_str());
                    if( !bOutputMode ) return false;
                }else{
                    m_EntrySpecified[optionIndex] = true;
                    int count = m_InputEntry[optionIndex].m_InputEntryValue.size();
                    if( m_DefinedEntry[optionIndex].m_InputOptionType == InputOptionType::MustAndOnce ||
                        m_DefinedEntry[optionIndex].m_InputOptionType == InputOptionType::CanBeLeftOut ){
                        if( count == 0 ){
                            text value = text::FromTextArray(&word, 1);
                            if( !AddInputValue(optionIndex, value, bPrint) ){
                                if( !bOutputMode ) return false;
                            }
                        }else{
                            if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) can be specified only once\n", word[0].c_str());
                            if( !bOutputMode ) return false;
                        }
                    }else if ( m_DefinedEntry[optionIndex].m_InputOptionType == InputOptionType::Multiple ||
                                m_DefinedEntry[optionIndex].m_InputOptionType == InputOptionType::Whatever ){
                        text value = text::FromTextArray(&word, 1);
                        if( !AddInputValue(optionIndex, value, bPrint) ){
                            if( !bOutputMode ) return false;
                        }
                    }
                }
            }
        }
        // check unspecified option value
        for(int i=0;i<m_DefinedEntry.size();i++){
            if( m_DefinedEntry[i].m_InputOptionType == InputOptionType::MustAndOnce ||
                m_DefinedEntry[i].m_InputOptionType == InputOptionType::Multiple ){
                if( m_InputEntry[i].m_InputEntryValue.size() == 0 ){
                    if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) must be specified at least once\n", m_DefinedEntry[i].m_Option.c_str());
                    if( !bOutputMode ) return false;
                }
            }
        }
        // add unspecified option value for CanBeLeftOut and Whatever
        for(int i=0;i<m_DefinedEntry.size();i++){
            if( m_DefinedEntry[i].m_InputOptionType == InputOptionType::CanBeLeftOut ||
                m_DefinedEntry[i].m_InputOptionType == InputOptionType::Whatever ){
                if( m_InputEntry[i].m_InputEntryValue.size() == 0 ){
                    AddInputValue(i, m_DefinedEntry[i].m_DefaultValue.m_ValueString, false);
                }
            }
        }
        return true;
    };
    bool InputParser::AddInputValue(int optionIndex, text value, bool bPrint){
        auto word = value.split();
        int requiredValue = m_DefinedEntry[optionIndex].m_NumberOfValues;
        if( m_DefinedEntry[optionIndex].m_InputNumberRequirement == InputNumberRequirement::Equal ){
            if( requiredValue != word.size() ){
                if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) requires exactly %d values but got %zd (%s)\n", 
                    m_DefinedEntry[optionIndex].m_Option.c_str(), requiredValue, word.size(), value.c_str());
                return false;
            }
        }else if( m_DefinedEntry[optionIndex].m_InputNumberRequirement == InputNumberRequirement::EqualOrGreaterThan ){
            if( requiredValue > word.size() ){
                if( bPrint ) fprintf(stderr, "Error (InputParser): (%s) requires no less than %d values but got %zd\n", 
                    m_DefinedEntry[optionIndex].m_Option.c_str(), requiredValue, word.size());
                return false;
            }
        }else if( m_DefinedEntry[optionIndex].m_InputNumberRequirement == InputNumberRequirement::Whatever ){
            // whatever
        }
        int count = m_InputEntry[optionIndex].m_InputEntryValue.size();
        m_InputEntry[optionIndex].m_InputEntryValue.resize(count + 1);
        m_InputEntry[optionIndex].m_InputEntryValue.back().m_ValueString = value;
        m_InputEntry[optionIndex].m_InputEntryValue.back().m_ValueVector = word;
        return true;
    }
};