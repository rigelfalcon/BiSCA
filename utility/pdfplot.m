function pdfplot(pd, x)
    %PDFPLOT Plots the PDF of input data
    %   data = input data vector
    %   title_str = title for the plot (optional)

    if nargin < 2
        xlim = pd.icdf([.001, .999]);
        x = linspace(xlim(1), xlim(2), 1001)';
    end

    pdf_x = pdf(pd, x);

    plot(x, pdf_x);

end