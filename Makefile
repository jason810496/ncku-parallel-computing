USER=jason810496
STUDENT=F74116720


init:
# copy $(HW)/TA2024_hw*.cpp to  $(HW)/$(STUDENT)_hw*.cpp
ifndef HW
	$(error HW is not set)
endif
# HW will be $-$, replace - with _
	touch $(HW)/$(STUDENT)_$(subst -,_,$(HW)).cpp
# replace mpicc with mpic++ in $(HW)/Makefile
	sed -i '' 's/mpicc/mpic++/g' $(HW)/Makefile
# replace TA2024 with $(STUDENT) in $(HW)/Makefile
	sed -i '' 's/TA2024_$(subst -,_,$(HW)).c/$(STUDENT)_$(subst -,_,$(HW)).cpp/g' $(HW)/Makefile
.PHONY: init

scp:
ifndef HW
	$(error HW is not set)
endif
# check if ~/$(HW) exists
	ssh ncku-parallel "test -d ~/$(HW) || mkdir ~/$(HW)"
ifndef FILE
	scp -r $(HW)/* ncku-parallel:~/$(HW)
else
	scp $(HW)/$(FILE) ncku-parallel:~/$(HW)
endif
.PHONY: scp

exe:
ifndef HW
	$(error HW is not set)
endif
ifndef CMD
	$(error CMD is not set)
endif
	ssh ncku-parallel 'cd ~/$(HW) && $(CMD)'
.PHONY: exe