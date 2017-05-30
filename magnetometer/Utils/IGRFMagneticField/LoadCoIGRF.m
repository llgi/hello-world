function [G,H] = LoadCoIGRF(year)
% function [G,H] = igrf_coeffs_gen(year);
% Enter year between 1900 - 2010 in multiple of 5 (e.g. 1915)
% From 2010-2015, enter year in multiple of 1 (e.g. 2011)

    maxYear = 2015;
    maxColumn = 26;

    load('g2015.txt');
    g = g2015; clear g2015
    load('h2015.txt');
    h = h2015; clear h2015

    G = zeros(14,14);
    H = zeros(14,14);

    year_to_column = [  1900	3
                        1905	4
                        1910	5
                        1915	6
                        1920	7
                        1925	8
                        1930	9
                        1935	10
                        1940	11
                        1945	12
                        1950	13
                        1955	14
                        1960	15
                        1965	16
                        1970	17
                        1975	18
                        1980	19
                        1985	20
                        1990	21
                        1995	22
                        2000	23
                        2005    24
                        2010    25
                        2015    26];

    if year > maxYear
        column = maxColumn;
    else
        indx = find(year_to_column >= year);
        column = year_to_column(indx(1),2);
    end

    svg = zeros(14,14); % Secular variation past year 2015
    svh = zeros(14,14); % Secular variation past year 2015

    if year > maxYear
        column = maxColumn;
        for i = 1:length(g)
            svg(g(i,1)+1,g(i,2)+1) = (year-maxYear)*g(i,26)/1e9;
        end
        for i = 1:length(h)
            svh(h(i,1)+1,h(i,2)+1) = (year-maxYear)*h(i,26)/1e9;
        end
    end

    for i = 1:length(g)
        G(g(i,1)+1,g(i,2)+1) = svg(g(i,1)+1,g(i,2)+1) + g(i,column)/1e9;
    end

    for i = 1:length(h)
        H(h(i,1)+1,h(i,2)+1) = svh(h(i,1)+1,h(i,2)+1) + h(i,column)/1e9;
    end
end