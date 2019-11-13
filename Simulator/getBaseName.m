function basename = getBaseName()
	stack = dbstack;
	basename = stack(2).name;
end