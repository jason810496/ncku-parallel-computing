CC := mpic++
CFLAGS := -Wall
TARGET := hw1_2
v := 1

all: $(TARGET)

$(TARGET):F74116720_hw1_2.cpp
	$(CC) -o $@ $^

judge: all
	@judge -v ${v} || printf "or \`make judge v=1\`"

clean:
	rm -f $(TARGET)
