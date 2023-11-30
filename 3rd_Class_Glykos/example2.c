#include <stdio.h>
#include <math.h>
#include <string.h>


int main()
{


int length;
int i;
float count = 0;
char seq [ 1000 ];


while (scanf("%s",seq) == 1 )
{
	length = strlen( seq );



	for (i=0 ; i < length ; i++)

{

	if ( seq[i] == 'G')
		{count = count + 1;}
 	if ( seq[i] == 'C')
		{count = count + 1;}
  

}


printf("The GC contect is %f\n",count);
count = 0;
}
}

