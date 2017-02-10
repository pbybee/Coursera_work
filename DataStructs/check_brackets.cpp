#include <iostream>
#include <stack>
#include <string>

struct Bracket {
    Bracket(char type, int position):
        type(type),
        position(position)
    {}

    bool Matchc(char c) {
        if (type == '[' && c == ']')
            return true;
        if (type == '{' && c == '}')
            return true;
        if (type == '(' && c == ')')
            return true;
        return false;
    }

    char type;
    int position;
};

int main() {

    std::string text;
    getline(std::cin, text);

    std::stack <Bracket> openBracket;
    int bracketFlag = 0;
    for (int position = 0; position < text.length(); ++position) {
        char next = text[position];
        Bracket newBracket(next,position);

        if (next == '(' || next == '[' || next == '{') {
            // Process opening bracket, write your code here
            openBracket.push(newBracket);
            bracketFlag = 1;
        }
        
//        Bracket *closeBracket;
        if (next == ')' || next == ']' || next == '}') {
            
            

            if (openBracket.empty() || !openBracket.top().Matchc(next)){
                std::cout << position+1;
                return(0);
            }
            if (!openBracket.empty() && openBracket.top().Matchc(next)){
                openBracket.pop();
            }
        } 
    }

    // Printing answer, write your code here
    if (openBracket.empty() && bracketFlag){
            std::cout << "Success";
    } else if (!openBracket.empty() && bracketFlag){
        std::cout << openBracket.top().position+1;
    }
    return 0;
}
