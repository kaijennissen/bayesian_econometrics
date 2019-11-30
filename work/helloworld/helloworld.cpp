#include <iostream>
#include <vector>
#include <string>
#include <experimental/optional>

using namespace std;

int main()
{
int a = 3;
int b = 5;
int c = a+b;
    vector<string> msg {"Hello", "C++", "World", "from", "VS Code!"};
    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
}