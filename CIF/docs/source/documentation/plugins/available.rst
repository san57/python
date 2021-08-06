#################
Available plugins
#################

List of all available plugins and requirements for each class:


	domain
		- CHIMERE, std
		- dummy, std
		- LMDZ, std


	chemistry
		- CHIMERE, gasJtab


	obsoperator
		- standard, std
			Requires:
				- model:
					- name:
					- version:
					- any: True
					- empty: False
				- obsvect:
					- name: standard
					- version: std
					- any: True
					- empty: True
				- controlvect:
					- name: standard
					- version: std
					- any: True
					- empty: True


	obsparser
		- WDCGG, std


	simulator
		- gausscost, std
			Requires:
				- controlvect:
					- name:
					- version:
					- any: True
					- empty: False
				- obsvect:
					- name:
					- version:
					- any: True
					- empty: False
				- obsoperator:
					- name: standard
					- version: std
					- any: True
					- empty: True
		- dummy_txt, std
			Requires:
				- controlvect:
					- name:
					- version:
					- any: True
					- empty: False
				- obsvect:
					- name:
					- version:
					- any: True
					- empty: False
				- obsoperator:
					- name: standard
					- version: std
					- any: True
					- empty: True


	controlvect
		- standard, std
			Requires:
				- domain:
					- name:
					- version:
					- any: True
					- empty: False
				- model:
					- name:
					- version:
					- any: True
					- empty: False


	measurements
		- random, std
			Requires:
				- domain:
					- name:
					- version:
					- any: True
					- empty: False
		- standard, std


	meteo
		- dummy, csv
			Requires:
				- domain:
					- name: dummy
					- version: std
					- any: False
					- empty: False
		- CHIMERE, std
			Requires:
				- domain:
					- name: CHIMERE
					- version: std
					- any: False
					- empty: False
		- LMDZ, mass-fluxes
			Requires:
				- domain:
					- name: LMDZ
					- version: std
					- any: False
					- empty: False


	platform
		- LSCE, obelix
			Requires:
				- model:
					- name:
					- version:
					- any: True
					- empty: True


	obsvect
		- standard, std
			Requires:
				- model:
					- name:
					- version:
					- any: True
					- empty: False
				- measurements:
					- name:
					- version:
					- any: True
					- empty: True


	mode
		- analytic, std
		- adj-tl_test, std
			Requires:
				- controlvect:
					- name: standard
					- version: std
					- any: True
					- empty: False
				- obsvect:
					- name: standard
					- version: std
					- any: True
					- empty: False
				- obsoperator:
					- name: standard
					- version: std
					- any: True
					- empty: True
		- 4dvar, std
			Requires:
				- controlvect:
					- name: standard
					- version: std
					- any: True
					- empty: False
				- obsvect:
					- name: standard
					- version: std
					- any: True
					- empty: False
				- obsoperator:
					- name: standard
					- version: std
					- any: True
					- empty: True
				- minimizer:
					- name: m1qn3
					- version: std
					- any: True
					- empty: True
				- simulator:
					- name: gausscost
					- version: std
					- any: True
					- empty: True
		- footprint, std
		- forward, std
			Requires:
				- controlvect:
					- name: standard
					- version: std
					- any: False
					- empty: True
				- obsoperator:
					- name: standard
					- version: std
					- any: False
					- empty: True
		- post-proc, std
			Requires:
				- controlvect:
					- name: standard
					- version: std
					- any: True
					- empty: True
				- obsvect:
					- name: standard
					- version: std
					- any: True
					- empty: True
				- obsoperator:
					- name: standard
					- version: std
					- any: True
					- empty: True


	model
		- dummy, std
			Requires:
				- meteo:
					- name: dummy
					- version: csv
					- any: False
					- empty: False
				- domain:
					- name: dummy
					- version: std
					- any: False
					- empty: False
				- fluxes:
					- name: dummy
					- version: nc
					- any: False
					- empty: True
		- LMDZ, std
			Requires:
				- meteo:
					- name: LMDZ
					- version: mass-fluxes
					- any: False
					- empty: False
				- domain:
					- name: LMDZ
					- version: std
					- any: False
					- empty: False
				- fluxes:
					- name: LMDZ
					- version: bin
					- any: False
					- empty: True
				- emis_species:
					- name: LMDZ
					- version: sflx
					- any: False
					- empty: True
				- chemistry:
					- name: CHIMERE
					- version: gasJtab
					- any: False
					- empty: False
		- CHIMERE, std
			Requires:
				- fluxes:
					- name: CHIMERE
					- version: AEMISSIONS
					- any: False
					- empty: True
				- domain:
					- name: CHIMERE
					- version: std
					- any: False
					- empty: False
				- chemistry:
					- name: CHIMERE
					- version: gasJtab
					- any: False
					- empty: False
				- meteo:
					- name: CHIMERE
					- version: std
					- any: False
					- empty: False


	fluxes
		- LMDZ, sflx
			Requires:
				- domain:
					- name: LMDZ
					- version: std
					- any: False
					- empty: False
		- dummy, txt
			Requires:
				- domain:
					- name: dummy
					- version: std
					- any: False
					- empty: False
		- dummy, nc
			Requires:
				- domain:
					- name: dummy
					- version: std
					- any: False
					- empty: False
		- LMDZ, bin
			Requires:
				- domain:
					- name: LMDZ
					- version: std
					- any: False
					- empty: False
		- CHIMERE, AEMISSIONS
			Requires:
				- domain:
					- name: CHIMERE
					- version: std
					- any: False
					- empty: False
				- chemistry:
					- name: CHIMERE
					- version: gasJtab
					- any: False
					- empty: False


	minimizer
		- congrad, std
			Requires:
				- simulator:
					- name: gausscost
					- version: std
					- any: True
					- empty: True
		- M1QN3, std
			Requires:
				- simulator:
					- name: gausscost
					- version: std
					- any: True
					- empty: True



