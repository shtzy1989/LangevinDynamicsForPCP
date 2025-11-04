#include "argumentParser.h"

ArgumentEntry::ArgumentEntry(bool positional, int numberOfValues, text optionShort, text optionLong, text helpText, std::vector<text> defaultValue, bool bRequired, ArgumentAction argumentAction){
    m_bPositional = positional;
    m_NumberOfValues = numberOfValues;
    m_OptionShort = optionShort;
    m_OptionLong = optionLong;
    m_HelpText = helpText;
    m_DefaultValue = defaultValue;
    m_Specified = false;
    m_ArgumentAction = argumentAction;
    m_bRequired = bRequired;
    m_Value.resize(1);
    if( !positional ){
        if( defaultValue.size() == 0 && numberOfValues != 0 ){
            fprintf(stderr, "Error: enough default values are not provided for optional argument\n");
            exit(0);
        }else if( defaultValue.size() == numberOfValues ){
            m_Value[0] = defaultValue;
        }else{
            fprintf(stderr, "Error: default values are provided but the number is different from expected value, so not used\n");
            exit(0);
        }
    }
}


void ArgumentParser::AddPositionalArgument(text optionShort, text helpText){
    m_PositionalEntryList.push_back(m_DefinedEntry.size());
    std::vector<text> empty;
    m_DefinedEntry.emplace_back(true, 1, optionShort, optionShort, helpText, empty, true, ArgumentAction::OnceOnly);
    // m_OptionalEntryMap[optionShort] = m_DefinedEntry.size();
    // m_OptionalEntryMap[optionShort] = m_DefinedEntry.size();
    m_TotalEntryMap[optionShort] = m_DefinedEntry.size();
    m_TotalEntryMap[optionShort] = m_DefinedEntry.size();
}
void ArgumentParser::AddOptionalArgument(text optionShort, text optionLong, text helpText, int numberofValues, std::vector<text> defaultValue, bool bRequired, ArgumentAction argumentAction){
    m_OptionalEntryList.push_back(m_DefinedEntry.size());
    m_DefinedEntry.emplace_back(false, numberofValues, optionShort, optionLong, helpText, defaultValue, bRequired, argumentAction);
    m_OptionalEntryMap[optionShort] = m_DefinedEntry.size();
    m_OptionalEntryMap[optionLong] = m_DefinedEntry.size();
    m_TotalEntryMap[optionShort] = m_DefinedEntry.size();
    m_TotalEntryMap[optionLong] = m_DefinedEntry.size();
}
bool ArgumentParser::Parse(int *argc, char **argv[]){
    if( !Parse0(argc, argv) ){
        PrintHelp();
        exit(0);
        return false;
    }else{
        return true;
    }
}
size_t ArgumentParser::GetNumberOfPositionalArgument(){
    return m_PositionalEntryList.size();
};
size_t ArgumentParser::GetNumberOfOptionalArgument(){
    return m_OptionalEntryList.size();
};
std::vector<text> ArgumentParser::GetOptionValue(text name){
    size_t index = m_TotalEntryMap[name] - 1;
    if( index != ULONG_MAX ){
        return m_DefinedEntry[index].m_Value[0];
    }else{
        fprintf(stderr, "Error: no argument \"%s\" defined\n", name.c_str());
        exit(0);
        std::vector<text> result;
        return result;
    }
};
std::vector<text> ArgumentParser::GetOptionValue(text name, int valueID){
    size_t index = m_TotalEntryMap[name] - 1;
    if( index != ULONG_MAX ){
        if( valueID >= 0 && valueID < m_DefinedEntry[index].m_Value.size() ){
            return m_DefinedEntry[index].m_Value[valueID];
        }else{
            fprintf(stderr, "Error: valueID is out of bound\n");
            exit(0);
            std::vector<text> result;
            return result;
        }
    }else{
        fprintf(stderr, "Error: no argument \"%s\" defined\n", name.c_str());
        exit(0);
        std::vector<text> result;
        return result;
    }
};
size_t ArgumentParser::GetNumberOfOptionValue(text name){
    size_t index = m_TotalEntryMap[name] - 1;
    if( index != ULONG_MAX ){
        if( m_DefinedEntry[index].m_Specified ){
            return m_DefinedEntry[index].m_Value.size();
        }else{
            return 0;
        }
    }else{
        fprintf(stderr, "Error: no argument \"%s\" defined\n", name.c_str());
        exit(0);
        return 0;
    }
};
bool ArgumentParser::GetOptionSpecified(text name){
    size_t index = m_TotalEntryMap[name] - 1;
    if( index != ULONG_MAX ){
        return m_DefinedEntry[index].m_Specified;
    }else{
        fprintf(stderr, "Error: no argument \"%s\" defined\n", name.c_str());
        exit(0);
        return false;
    }
};
void ArgumentParser::PrintAllValues(){
    for(int i=0;i<m_PositionalEntryList.size();i++){
        int argumentIndex = m_PositionalEntryList[i];
        fprintf(stderr, "Positional %s %3d %-10s ", m_DefinedEntry[argumentIndex].m_Specified ? "T" : "F", i, m_DefinedEntry[argumentIndex].m_OptionShort.c_str());
        for(int k=0;k<m_DefinedEntry[argumentIndex].m_Value.size();k++){
            for(int j=0;j<m_DefinedEntry[argumentIndex].m_Value[k].size();j++){
                fprintf(stderr, "%-10s ", m_DefinedEntry[argumentIndex].m_Value[k][j].c_str());
            }
        }
        fprintf(stderr, "\n");
    }
    for(int i=0;i<m_OptionalEntryList.size();i++){
        int argumentIndex = m_OptionalEntryList[i];
        fprintf(stderr, "Optional   %s %3d %-10s ", m_DefinedEntry[argumentIndex].m_Specified ? "T" : "F", i, m_DefinedEntry[argumentIndex].m_OptionShort.c_str());
        for(int k=0;k<m_DefinedEntry[argumentIndex].m_Value.size();k++){
            if( k != 0 ){
                fprintf(stderr, ": ");
            }
            for(int j=0;j<m_DefinedEntry[argumentIndex].m_Value[k].size();j++){
                fprintf(stderr, "%-10s ", m_DefinedEntry[argumentIndex].m_Value[k][j].c_str());
            }
        }
        fprintf(stderr, "\n");
    }
}
void ArgumentParser::PrintHelp(){
    fprintf(stderr, "command usage\n");
    for(int i=0;i<m_PositionalEntryList.size();i++){
        int argumentIndex = m_PositionalEntryList[i];
        fprintf(stderr, "    Positional %-20s            %-s\n", 
            m_DefinedEntry[argumentIndex].m_OptionShort.c_str(),
            m_DefinedEntry[argumentIndex].m_HelpText.c_str());
    }
    for(int i=0;i<m_OptionalEntryList.size();i++){
        int argumentIndex = m_OptionalEntryList[i];
        fprintf(stderr, "    Optional   %-5s %-14s %8s %1d %-20s Default: ", 
            m_DefinedEntry[argumentIndex].m_OptionShort.c_str(),
            m_DefinedEntry[argumentIndex].m_OptionLong.c_str(),    
            m_DefinedEntry[argumentIndex].m_bRequired ? "Required" : "",            
            m_DefinedEntry[argumentIndex].m_NumberOfValues,
            m_DefinedEntry[argumentIndex].m_HelpText.c_str());
        for(int j=0;j<m_DefinedEntry[argumentIndex].m_DefaultValue.size();j++){
            fprintf(stderr, "%10s ", m_DefinedEntry[argumentIndex].m_DefaultValue[j].c_str());
        }
        fprintf(stderr, "\n");
    }
}
bool ArgumentParser::Parse0(int *argc, char **argv[]){
    std::vector<text> argumentList;
    if( *argc > 0 ){
        argumentList.reserve(*argc - 1);
    }
    for(int i=1;i<*argc;i++){
        argumentList.push_back((*argv)[i]);
    }
    if( argumentList.size() == 1 && ( argumentList[0].ContentEqual("-h") || argumentList[0].ContentEqual("--help") ) ){
        PrintHelp();
        exit(0);
    }
    for(int i=0;i<argumentList.size();i++){
        text argument = argumentList[i];
        if( i < m_PositionalEntryList.size() ){
            // positional arguments
#if __cplusplus >= 202002L
            if( m_OptionalEntryMap.contains(argument) ){
#else
            size_t index = m_OptionalEntryMap[argument] - 1;
            if( index != ULONG_MAX ){
#endif
                fprintf(stderr, "Error: positional argument is expected but \"%s\" is specified\n", argument.c_str());
                return false;
            }else{
                int argumentIndex = m_PositionalEntryList[i];
                m_DefinedEntry[argumentIndex].m_Value[0].push_back(argument);
                m_DefinedEntry[argumentIndex].m_Specified = true;
            }
        }else{
            // optional arguments
            size_t argumentIndex = m_OptionalEntryMap[argument] - 1;
            if( argumentIndex != ULONG_MAX ){
                // defined argument
                if( m_DefinedEntry[argumentIndex].m_Specified ){
                    // from the second time this argument is defined
                    if( m_DefinedEntry[argumentIndex].m_ArgumentAction == ArgumentAction::OnceOnly ){
                        fprintf(stderr, "Error: optional argument \"%s\" can only be specified once\n", 
                            m_DefinedEntry[argumentIndex].m_OptionShort.c_str());
                        return false;
                    }else if( m_DefinedEntry[argumentIndex].m_ArgumentAction == ArgumentAction::Overwrite ){
                        // m_DefinedEntry[argumentIndex].m_Value[0].clear();
                        // nothing
                    }else if( m_DefinedEntry[argumentIndex].m_ArgumentAction == ArgumentAction::Push ){
                        m_DefinedEntry[argumentIndex].m_Value.resize(m_DefinedEntry[argumentIndex].m_Value.size() + 1);
                        m_DefinedEntry[argumentIndex].m_Value.back().resize(m_DefinedEntry[argumentIndex].m_NumberOfValues);
                    }
                }
                bool bOK = true;
                if( i + m_DefinedEntry[argumentIndex].m_NumberOfValues < argumentList.size() ){
                    for(int j=0;j<m_DefinedEntry[argumentIndex].m_NumberOfValues;j++){
#if __cplusplus >= 202002L      
                        if( m_OptionalEntryMap.contains(argumentList[i+j+1]) ){
#else
                        size_t index = m_OptionalEntryMap[argumentList[i+j+1]] - 1;
                        if( index != ULONG_MAX ){
#endif
                            // premature option list for this argument
                            bOK = false;
                            break;
                        }
                        m_DefinedEntry[argumentIndex].m_Value.back()[j] = argumentList[i+j+1];
                    }
                    m_DefinedEntry[argumentIndex].m_Specified = true;
                    i = i + m_DefinedEntry[argumentIndex].m_NumberOfValues;
                }else{
                    bOK = false;
                }
                if( !bOK ){
                    fprintf(stderr, "Error: optional argument \"%s\" expects %d argument(s)\n", 
                        m_DefinedEntry[argumentIndex].m_OptionShort.c_str(), 
                        m_DefinedEntry[argumentIndex].m_NumberOfValues);
                    return false;
                }
            }else{
                // undefined argument
                fprintf(stderr, "Error: undefined option %s\n", argument.c_str());
                return false;
            }
        }
    }
    for(int i=0;i<m_PositionalEntryList.size();i++){
        int argumentIndex = m_PositionalEntryList[i];
        if( m_DefinedEntry[argumentIndex].m_Value[0].size() != 1 ){
            fprintf(stderr, "Error: expecting %zd positional argument but only %d specified\n", 
                m_PositionalEntryList.size(),
                i);
            return false;
        }
    }
    for(int i=0;i<m_OptionalEntryList.size();i++){
        int argumentIndex = m_OptionalEntryList[i];
        if( m_DefinedEntry[argumentIndex].m_bRequired && !m_DefinedEntry[argumentIndex].m_Specified ){
            fprintf(stderr, "Error: optional argument %s is marked as required but not specified\n", 
                m_DefinedEntry[argumentIndex].m_OptionShort.c_str());
            return false;
        }
    }
    
    return true;
};
