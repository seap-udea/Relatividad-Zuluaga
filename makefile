include mybook.conf

test:
	@echo "Testing installation..."	
	@jupyter nbconvert --to notebook --execute nb/$(TESTS)

clean:
	@echo "Cleaning crap..."
	@-find . -name "*~" -exec rm -rf {} \;
	@-find . -name "#*#" -delete
	@-find . -name "#*" -delete
	@-find . -name ".#*" -delete
	@-find . -name ".#*#" -delete
	@-find . -name "*.nbconvert.*" -delete
	@-find . -name ".DS_Store" -delete
	@-find . -name ".Icon*" -delete
	@-find . -name '__pycache__' -type d | xargs rm -fr

cleanall:clean
	@echo "Cleaning all..."
	@-rm -rf nb/* html/*

update:
	@echo "Adding files..."
	@-git add -Af

commit:clean update
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin master

pull:
	@echo "Pulling new files..."
	@-git reset --hard HEAD
	@-git pull origin master

pack:
	@echo "Packing data..."
	@bash .store/pack.sh pack

unpack:
	@echo "Unpacking data..."
	@bash .store/pack.sh unpack

