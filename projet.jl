using IntervalArithmetic, IntervalRootFinding


#----------------------------FONCTION F ET SA DÉRIVÉE-----------------------------
#fonction f étudiée
function f(x)
    return x^3 + 2*x^2 - 5*x -6
end

#dérivée de f
function df(x)
    return 3*x^2 + 4*x - 5
end
#---------------------------------------------------------------------------------



#--------------------------- QUESTION 2 -------------------------------------------
function F(X::Interval)
    if(!isempty(X))
        x1 = (-4-sqrt(76))/6
        x2 = (-4+sqrt(76))/6
        
        if X.hi <= x1
            return Interval(f(X.lo), f(X.hi))
        else
            if X.hi <= x2 
                if X.lo <= x1
                    return Interval(min(f(X.lo),f(X.hi)), f(x1))
                else
                    return Interval(f(X.hi), f(X.lo))
                end
            else #X.hi > x2
                if X.lo <= x1
                    return Interval(min(f(X.lo),f(x2)), max(f(X.hi),f(x1)))
                else
                    if X.lo <= x2
                        return Interval(f(x2), max(f(X.lo),f(X.hi)))
                    else
                        return Interval(f(X.lo), f(X.hi))
                    end
                end
            end
        end
    else
        return X
    end
end
#---------------------------------------------------------------------------------



#------------------------------- QUESTION 3 --------------------------------------
function dF(X::Interval)
    if(!isempty(X))
        x0 = -2/3
        if X.hi <= x0
            return Interval(df(X.hi), df(X.lo))
        else
        	if X.lo <= x0
        		return Interval(df(x0), max(df(X.lo),df(X.hi)))
        	else #X.lo >= x0
        		return Interval(df(X.lo), df(X.hi)) 
        	end
        end
    else
        return X   
    end
end
#---------------------------------------------------------------------------------



#-------------------------------- QUESTION 4 -------------------------------------
function Q4(X::Interval, Y::Interval)
    resInf = X.lo
    resSup = X.hi
	#résolution de inf(Y) = f(x) et sup(Y) = f(x) dans X
	println("debut newton")
	solYinf = roots(x -> x^3 + 2*x^2 - 5*x -6 - Y.lo, X)
	solYsup = roots(x -> x^3 + 2*x^2 - 5*x -6 - Y.hi, X)
	println("fin newton")
	solutions = append!([i.interval for i in solYinf], [i.interval for i in solYsup])
	
	#s'il y a des solutions dans X, on regarder si on peut contracter X
	if !isempty(solutions) 
	    #on essaye de contracter sur inf(X)
	    argSolXinf = argmin([i.lo for i in solutions])
	    solXinf = solutions[argSolXinf]
	    if argSolXinf <= length(solYinf) #solution de f(x) = inf(Y)
	        #si f croit en solXinf : on peut contracter
	        dFx = dF(solXinf)
	        if dFx.lo > 0 
	            resInf = solXinf.lo
	        end
	    else #solution de f(x) = sup(Y)
	        #si f décroit en solXinf : on peut contracter
	        dFx = dF(solXinf)
	        if dFx.hi < 0 
	            resInf = solXinf.lo
	        end
	    end
	    #on essaye de contracter sur sup(X)
	    argSolXsup = argmax([i.hi for i in solutions])
	    solXsup = solutions[argSolXsup]
	    if argSolXsup <= length(solYsup) #solution de f(x) = inf(Y)
	        #si f décroit en solXsup : on peut contracter
	        dFx = dF(solXsup)
	        if dFx.sup < 0 
	            resSup = solXsup.hi
	        end
	    else #solution de f(x) = sup(Y)
	        #si f croit en solXsup : on peut contracter
	        dFx = dF(solXsup)
	        if dFx.lo > 0 
	            resSup = solXsup.hi
	        end
	    end
	    
	else #pas de solutions dans X
	    #si F(X) n'est pas inclus dans Y on retourne l'intervalle vide
	    if !((F(X).lo >= Y.lo) & (F(X).hi <= Y.hi))
	        return emptyinterval()
	    end #sinon on retourne X
	end
	
	return Interval(resInf,resSup)
end
#---------------------------------------------------------------------------------



#------------------------------ QUESTION 6 ---------------------------------------
function Q6(X::Interval, Y::Interval)
    #alterner contraction sur X (Q4) et sur Y (Q2)
    hullConsistent = false
    while !hullConsistent #on contracte tant qu'on est pas hull consistent
        #contraction sur Y
        Y2 = intersect(F(X), Y)
        
        #contraction sur X
        X2 = intersect(Q4(X,Y), X)
        
        #on vérifie si on a contracté
        if (Y2 == Y) & (X2 == X)
            hullConsistent = true
        end
        X = X2
        Y = Y2
    end
    
    return IntervalBox(X,Y)
end
#---------------------------------------------------------------------------------



#----------------------------- NEWTON --------------------------------------------

#---------------------------------------------------------------------------------



#------------------------------- TESTS -------------------------------------------

function main()
	X1 = Interval(-6,-3.5)  
	X2 = Interval(-3,5)  
	X3 = Interval(3,5)   
	X4 = Interval(0,0)   
	X5 = emptyinterval() 
	X6 = Interval(-2,2)
	X7 = Interval(0,3)
	X8 = Interval(-3,3)         
	
	println("tests extension F(x) : ")
	println(F(X1))
	println(F(X2))
	println(F(X3))
	println(F(X4))
    println(F(X5))
    println(F(X6))
    println(F(X7))
    println(F(X8))
    println()
	
	println("tests extension dF(x) : ")
	println(dF(X1))
	println(dF(X2))
	println(dF(X3))
	println(dF(X4))
	println(dF(X5))
	println(dF(X6))
    println(dF(X7))
    println(dF(X8))
    println()
	
	Y1 = Interval(-2,2)
	Y2 = Interval(-8,-3)
	Y3 = Interval(0,9)
	Y4 = Interval(-9,10)
	Y5 = emptyinterval()
	Y6 = Interval(-9,6)
	Y7 = Interval(-6,6)
	Y8 = Interval(-6, f((-4-sqrt(76))/6))
	
	println("tests question 4 : ")	
	println(Q4(X1,Y1))
	println(Q4(X2,Y2))
	println(Q4(X7,Y3))
	println(Q4(X4,Y4))
	println(Q4(X4,Y5))
	println(Q4(X5,Y2))
	println(Q4(X6,Y6))
	println(Q4(X7,Y7))
	println(Q4(X8,Y8))
	println()
	
	println("tests question 6 : ")
	println(Q6(X1,Y1))
	println(Q6(X2,Y2))
	println(Q6(X3,Y3))
	println(Q6(X4,Y4))
	println(Q6(X4,Y5))
	println(Q6(X5,Y2))
	println(Q6(X6,Y6))
	println(Q6(X7,Y7))
	println(Q6(X8,Y8))
	println()
end

main()

