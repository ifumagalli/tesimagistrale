tesimagistrale
==============

Codici Fenics-Python per la tesi e il progetto di PACS

==============

bin\Lab01\, bin\Lab02\: codici dei laboratori di CSFluid 2013/14 scritti dal Dott. Mattia Penati.

bin\demo_optimization.py: dalle demo di dolfin, per avere ispirazione.

bin\tests\: test per vedere perché non funzionava la valutazione di un Function in un punto sulle macchine del Dipartimento: nomefile.out è l'output di $> python nomefile.py

bin\myschur.py[_old]: versioni di bin\Lab02\schur.py che vengono utilizzate; l'_old è come schur.py e necessita di PETSc, il non _old funziona ovunque, anche se è meno ottimizzato

==============

Branches:
	sincro: contiene solo due commit in più, che avevo fatto per sincronizzare ad un ricevimento, ma non servono più (da eliminare tra qualche giorno)

==============

Struttura descrizioni commit (ogni sezione può esserci o meno):
	<descrizione delle modifiche importanti>
	TODO:
	<cose da fare nel commit successivo>
	DONE:
	<conferma di aver evaso i TODO del commit precedente>
	MINOR:
	<modifiche non sostanziali>
	Note:
	<varie ed eventuali>
"Modifiche trascurabili in <file>" = <file> presenta solo modifiche trascurabili o modifiche dipendenti da modifiche sostanziali segnalate precedentemente.
Le modifiche di ogni file, se non trascurabili, devono essere tutte segnalate, altrimenti si metta il file in Modifiche trascurabili.
