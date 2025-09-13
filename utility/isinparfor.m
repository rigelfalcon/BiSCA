function ip = isinparfor()
    job = getCurrentJob();
    ip = ~isempty(job);
end