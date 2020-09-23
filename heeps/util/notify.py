import os

def notify(send_message, send_to, send_subject='HEEPS noreply'):

    ''' send a message by email. Subject defaults to "HEEPS noreply". '''

    if send_to is not None:
        os.system('echo "%s" | mail -s "%s" %s'%(send_message, send_subject, send_to))
