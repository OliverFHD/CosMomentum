#include "assert.h"
using namespace std;

namespace error_handling{

enum VALUE_TEST_MODUS {EQUAL, NOT_EQUAL, SMALLER, LARGER, SMALLER_EQUAL, LARGER_EQUAL};


void general_error_message(string error_message){
  cerr << error_message << '\n';
  exit(1);
}

void general_warning(string warning_message){
  cout << warning_message << '\n';
}


/*
 * Tests whether value x fullfils a certain criterion, defined by the value x_criterion and the nature of the criterion specified by MODUS.
 * If criterion is not met, an error message is displayed and the program abborts.
 * 
 */
 void test_value_of_double(double x, double x_criterion, VALUE_TEST_MODUS MODUS, string error_message){
   
   switch(MODUS){
     case EQUAL:
       if(x != x_criterion){
         cerr << error_message << '\n';
         assert(x = x_criterion);
       }
       break;
     case NOT_EQUAL:
       if(x == x_criterion){
         cerr << error_message << '\n';
         assert(x != x_criterion);
       }
       break;
     case SMALLER:
       if(x >= x_criterion){
         cerr << error_message << '\n';
         assert(x < x_criterion);
       }
       break;
     case LARGER:
       if(x <= x_criterion){
         cerr << error_message << '\n';
         assert(x > x_criterion);
       }
       break;
     case SMALLER_EQUAL:
       if(x > x_criterion){
         cerr << error_message << '\n';
         assert(x <= x_criterion);
       }
       break;
     case LARGER_EQUAL:
       if(x < x_criterion){
         cerr << error_message << '\n';
         assert(x >= x_criterion);
       }
       break;
   }
   
 }

/*
 * Tests whether value x fullfils a certain criterion, defined by the value x_criterion and the nature of the criterion specified by MODUS.
 * If criterion is not met, an error message is displayed and the program abborts.
 * 
 */
 void test_value_of_int(int x, int x_criterion, VALUE_TEST_MODUS MODUS, string error_message){
   
   switch(MODUS){
     case EQUAL:
       if(x != x_criterion){
         cerr << error_message << '\n';
         assert(x = x_criterion);
       }
       break;
     case NOT_EQUAL:
       if(x == x_criterion){
         cerr << error_message << '\n';
         assert(x != x_criterion);
       }
       break;
     case SMALLER:
       if(x >= x_criterion){
         cerr << error_message << '\n';
         assert(x < x_criterion);
       }
       break;
     case LARGER:
       if(x <= x_criterion){
         cerr << error_message << '\n';
         assert(x > x_criterion);
       }
       break;
     case SMALLER_EQUAL:
       if(x > x_criterion){
         cerr << error_message << '\n';
         assert(x <= x_criterion);
       }
       break;
     case LARGER_EQUAL:
       if(x < x_criterion){
         cerr << error_message << '\n';
         assert(x >= x_criterion);
       }
       break;
   }
   
 }


/*
 * Tests whether value x fullfils a certain criterion, defined by the value x_criterion and the nature of the criterion specified by MODUS.
 * If criterion is not met, a warning message is displayed.
 * 
 */
 void warning_value_of_double(double x, double x_criterion, VALUE_TEST_MODUS MODUS, string warning_message){
   
   switch(MODUS){
     case EQUAL:
       if(x != x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case NOT_EQUAL:
       if(x == x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case SMALLER:
       if(x >= x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case LARGER:
       if(x <= x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case SMALLER_EQUAL:
       if(x > x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case LARGER_EQUAL:
       if(x < x_criterion){
         cout << warning_message << '\n';
       }
       break;
   }
   
 }


/*
 * Tests whether value x fullfils a certain criterion, defined by the value x_criterion and the nature of the criterion specified by MODUS.
 * If criterion is not met, a warning message is displayed.
 * 
 */
 void warning_value_of_int(int x, int x_criterion, VALUE_TEST_MODUS MODUS, string warning_message){
   
   switch(MODUS){
     case EQUAL:
       if(x != x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case NOT_EQUAL:
       if(x == x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case SMALLER:
       if(x >= x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case LARGER:
       if(x <= x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case SMALLER_EQUAL:
       if(x > x_criterion){
         cout << warning_message << '\n';
       }
       break;
     case LARGER_EQUAL:
       if(x < x_criterion){
         cout << warning_message << '\n';
       }
       break;
   }
   
 }

}
