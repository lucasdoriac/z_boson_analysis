#include <iostream>
#include <bitset>
using namespace std;

//My triggers
enum Trigger{

    JUMP   = 1 << 0,
    ATTACK = 1 << 1,
    DASH   = 1 << 2,
    SHIELD = 1 << 3
};

void aula_1();
void aula_2();

int main()
{
    aula_1();
}

void aula_1(){
    
    int triggers = 0;

    triggers |= JUMP;
    triggers |= ATTACK;
    triggers |= SHIELD;

    cout << "Value of triggers = " << triggers << endl;

    //We got 11. Is this reasonable? Its 2^0 + 2^1 + 0 + 2^3 = 1 + 2 + 0 + 8 = 11.
    //But what if i want to print its bin form?
    //Podemos usar o std::bitset<n>, onde n = número da representação de bits:
    cout << "Value of triggers in bits: " << std::bitset<4>(triggers) << endl;

    triggers ^= ATTACK;
    triggers &= ~JUMP;

    //So triggers now will have the value 1000 = 8.

    cout << "Value of triggers = " << triggers << endl;

}

void aula_2(){
    //Construindo um sistema de triggers.

    

}