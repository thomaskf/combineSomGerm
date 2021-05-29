all : combineSomGerm

combineSomGerm : combineSomGerm.cpp combineSomGerm.h
	g++ -o combineSomGerm combineSomGerm.cpp

clean :
	rm combineSomGerm
