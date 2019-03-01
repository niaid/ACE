# Locus

## Login via ssh

```bash
ssh -X username@ai-submit2.niaid.nih.gov
## enter password
```



## Interactive session

```bash
qrsh -l h_vmem=32G
```



## Access folders via sftp

- Cyberduck -> Open Connection (icon top left of the window)
![](assets/img/cyberduck.png)





## Copy folders from Locus to Laptop via command line

```bash
## example command: replace username with your username!
scp -r username@ai-submit2.niaid.nih.gov:/classhome/username ~/Desktop
```



## More information for the curious

[Locus website](https://locus.niaid.nih.gov) - if you don't have a Locus server account, the first time you log into the website (using your regular NIH creds), you will get an email telling you how to request an account.

[Request an account](https://locus.niaid.nih.gov/userportal/documentation.php#Getting-Started/Request-an-Account) - in case you lose the email

- If you are curious how we came up with the scp command:
```bash
 ## basic command structure
 ## `-r` means "recursive" so we copy the folder and everything in it.
 scp -r fromfoldername tofoldername

 ## general command
 scp -r username@servername:serverfolderpath laptopfolderpath
```

