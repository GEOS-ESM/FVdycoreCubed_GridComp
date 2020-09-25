void putenv_(str, strl) 
char *str; 
long strl; 
{ 
        void *malloc(); 
        char *putstr = malloc(strl+1); 
        strncpy(putstr, str, strl); 
        putstr[strl] = 0; 
        putenv(putstr); 
}
